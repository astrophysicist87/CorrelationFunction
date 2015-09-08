#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

#include "CFWR.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

void CorrelationFunction::Get_GF_HBTradii(FO_surf* FOsurf_ptr, int folderindex)
{
   bool includezeroes = false;

   for(int iKT = 0; iKT < n_localp_T; iKT++)
   {
      //cout << "\t Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      *global_out_stream_ptr << "   - Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
      for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
      {
         *global_out_stream_ptr << "\t\t --> Calculating K_phi = " << K_phi[iKphi] << " ..." << endl;
         Reset_EmissionData();
         SetEmissionData(FOsurf_ptr, K_T[iKT], K_phi[iKphi], includezeroes);
	 if (corrfuncdim == 1)
	 {
		Cal_correlationfunction_1D();
		Fit_Correlationfunction1D('o', iKT, iKphi);
		Fit_Correlationfunction1D('s', iKT, iKphi);
		Fit_Correlationfunction1D('l', iKT, iKphi);
		Output_Correlationfunction_1D(iKT, iKphi, folderindex);
	 }
	 else if (corrfuncdim == 3)
	 {
		//Cal_correlationfunction_3D(iKT, iKphi);
		//if (lambdaflag)
		//	Fit_Correlationfunction3D_withlambda(iKT, iKphi);
		//else Fit_Correlationfunction3D(iKT, iKphi);
		//Output_Correlationfunction_3D(iKT, iKphi, folderindex);
        cerr << "corrfuncdim == 3 not currently supporter!" << endl;
	 }
      }
      R2_Fourier_transform(iKT, global_plane_psi);
   }
}

void CorrelationFunction::Cal_correlationfunction_1D()
{
    // Get spectra first
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
    {
        double integ1 = dN_dypTdpTdphi_moments[target_particle_id][0][0][0][ipt][ipphi];
        double integ2 = dN_dypTdpTdphi_moments[target_particle_id][0][0][1][ipt][ipphi];
        spectra[target_particle_id][ipt][ipphi] = integ1*integ1 + integ2*integ2;
    }

    // Then compute full correlation function
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
    for (int iq = 0; iq < qnpts; ++iq)
    for (int iqax = 0; iqax < 3; ++iqax)
    {
        double integ1 = dN_dypTdpTdphi_moments[target_particle_id][iq][iqax][0][ipt][ipphi];
        double integ2 = dN_dypTdpTdphi_moments[target_particle_id][iq][iqax][1][ipt][ipphi];
        CFvals[ipt][ipphi][iq][iqax] = (integ1*integ1 + integ2*integ2) / spectra[target_particle_id][ipt][ipphi];
    }

    return;
}

int CorrelationFunction::Read_correlationfunction_1D(int iKT, int iKphi)
{
   if(fabs(K_y) > 1e-16)
   {
       //cout<<"not support for y not equals 0 yet!" << endl;
       return (0);
   }

   double local_K_T = K_T[iKT];
   double localK_phi = K_phi[iKphi];

ostringstream corrfn_stream;
corrfn_stream << global_path;

if (fabs(local_K_T) < 1e-16) corrfn_stream << "/correlfunct1D" << "_" << particle_name << "_kt_0.0_phi_";
else if (fabs(local_K_T - 1.) < 1e-16) corrfn_stream << "/correlfunct1D" << "_" << particle_name << "_kt_1.0_phi_";
else corrfn_stream << "/correlfunct1D" << "_" << particle_name << "_kt_" << local_K_T << "_phi_";

if (fabs(localK_phi) < 1e-16) corrfn_stream << "0.0.dat";
else corrfn_stream << localK_phi << ".dat";

if ( !fexists( corrfn_stream.str().c_str() ) )
{
	cerr << corrfn_stream.str().c_str() << ": File not found" << endl;
	*global_out_stream_ptr << corrfn_stream.str().c_str() << ": File not found" << endl;
	return (1);
}
  ifstream corrfn(corrfn_stream.str().c_str());

   for(int i = 0; i < qnpts; i++)
   {
            corrfn >> q_out[i];
            corrfn >> Correl_1D_out[i];
            corrfn >> Correl_1D_out_err[i];
            corrfn >> q_side[i];
            corrfn >> Correl_1D_side[i];
            corrfn >> Correl_1D_side_err[i];
            corrfn >> q_long[i];
            corrfn >> Correl_1D_long[i];
            corrfn >> Correl_1D_long_err[i];
   }

   return (0);
}

//*********************************************************************
// Functions used for multidimension fit
void CorrelationFunction::Fit_Correlationfunction1D(char osl_switch, int iKT, int iKphi)
{
  const int data_length = qnpts;  // # of points
  const size_t n_para = 2;  // # of parameters

  // allocate space for a covariance matrix of size p by p
  gsl_matrix *covariance_ptr = gsl_matrix_alloc (n_para, n_para);

  // allocate and setup for generating gaussian distibuted random numbers
  gsl_rng_env_setup ();
  const gsl_rng_type *type = gsl_rng_default;
  gsl_rng *rng_ptr = gsl_rng_alloc (type);

  //set up test data
  struct Correlationfunction1D_data Correlfun1D_data;
  Correlfun1D_data.data_length = data_length;
  Correlfun1D_data.q = new double [data_length];
  Correlfun1D_data.y = new double [data_length];
  Correlfun1D_data.sigma = new double [data_length];

switch(osl_switch)
{
  case 'o':
  //cout << "in the out loop!" << endl;
  for(int i=0; i<data_length; i++)
  {
     Correlfun1D_data.q[i] = q_out[i];
     Correlfun1D_data.y[i] = Correl_1D_out[i];
     Correlfun1D_data.sigma[i] = Correl_1D_out_err[i];
  }
  break;
  case 's':
  //cout << "in the side loop!" << endl;
  for(int i=0; i<data_length; i++)
  {
     Correlfun1D_data.q[i] = q_side[i];
     Correlfun1D_data.y[i] = Correl_1D_side[i];
     Correlfun1D_data.sigma[i] = Correl_1D_side_err[i];
  }
  break;
  case 'l':
  //cout << "in the long loop!" << endl;
  for(int i=0; i<data_length; i++)
  {
     Correlfun1D_data.q[i] = q_long[i];
     Correlfun1D_data.y[i] = Correl_1D_long[i];
     Correlfun1D_data.sigma[i] = Correl_1D_long_err[i];
  }
  break;
}

  double para_init[n_para] = {1.0, 1.0};  // initial guesse of parameters

  gsl_vector_view xvec_ptr = gsl_vector_view_array (para_init, n_para);
  
  // set up the function to be fit 
  gsl_multifit_function_fdf target_func;
  target_func.f = &Fittarget_correlfun1D_f;        // the function of residuals
  target_func.df = &Fittarget_correlfun1D_df;      // the gradient of this function
  target_func.fdf = &Fittarget_correlfun1D_fdf;    // combined function and gradient
  target_func.n = data_length;              // number of points in the data set
  target_func.p = n_para;              // number of parameters in the fit function
  target_func.params = &Correlfun1D_data;  // structure with the data and error bars

  const gsl_multifit_fdfsolver_type *type_ptr = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver *solver_ptr 
       = gsl_multifit_fdfsolver_alloc (type_ptr, data_length, n_para);
  gsl_multifit_fdfsolver_set (solver_ptr, &target_func, &xvec_ptr.vector);

  size_t iteration = 0;         // initialize iteration counter
  //print_fit_state_1D (iteration, solver_ptr);
  int status;  		// return value from gsl function calls (e.g., error)
  do
  {
      iteration++;
      
      // perform a single iteration of the fitting routine
      status = gsl_multifit_fdfsolver_iterate (solver_ptr);

      // print out the status of the fit
      //cout << "status = " << gsl_strerror (status) << endl;

      // customized routine to print out current parameters
      //print_fit_state_1D (iteration, solver_ptr);

      if (status)    // check for a nonzero status code
      {
          break;  // this should only happen if an error code is returned 
      }

      // test for convergence with an absolute and relative error (see manual)
      status = gsl_multifit_test_delta (solver_ptr->dx, solver_ptr->x, 
                                        fit_tolarence, fit_tolarence);
  }
  while (status == GSL_CONTINUE && iteration < fit_max_iterations);

  // calculate the covariance matrix of the best-fit parameters
  gsl_multifit_covar (solver_ptr->J, 0.0, covariance_ptr);

  // print out the covariance matrix using the gsl function (not elegant!)
  //cout << endl << "Covariance matrix: " << endl;
  //gsl_matrix_fprintf (stdout, covariance_ptr, "%g");

  //cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
  //cout.precision (5);		                // # of digits in doubles

  int width = 7;		// setw width for output
  //cout << endl << "Best fit results:" << endl;
  //cout << "lambda = " << setw (width) << get_fit_results (0, solver_ptr)
  //  << " +/- " << setw (width) << get_fit_err (0, covariance_ptr) << endl;

  //cout << "R2 = " << setw (width) << get_fit_results (1, solver_ptr)
  //  << " +/- " << setw (width) << get_fit_err (1, covariance_ptr) << endl;

  //cout << "status = " << gsl_strerror (status) << endl;
  //cout << "------------------------------------------------------------------" << endl;

  double chi = gsl_blas_dnrm2(solver_ptr->f);
  double dof = data_length - n_para;
  double c = GSL_MAX_DBL(1, chi/sqrt(dof));

  lambda_Correl[iKT][iKphi] = get_fit_results(0, solver_ptr);
  lambda_Correl_err[iKT][iKphi] = c*get_fit_err(0, covariance_ptr);
  //cout << "final results: " << endl;
  //cout << scientific << setw(10) << setprecision(5) 
  //     << "chisq/dof = " << chi*chi/dof << endl;
  //cout << scientific << setw(10) << setprecision(5) 
  //     << " lambda = " << lambda_Correl << " +/- " << lambda_Correl_err << endl;
switch(osl_switch)
{
	case 'o':
		R2_out[iKT][iKphi] = get_fit_results(1, solver_ptr)*hbarC*hbarC;
		R2_out_err[iKT][iKphi] = c*get_fit_err(1, covariance_ptr)*hbarC*hbarC;
		//cout << " R2_out = " << R2_out_Correl << " +/- " << R2_out_Correl_err << endl;
		break;
	case 's':
		R2_side[iKT][iKphi] = get_fit_results(1, solver_ptr)*hbarC*hbarC;
		R2_side_err[iKT][iKphi] = c*get_fit_err(1, covariance_ptr)*hbarC*hbarC;
		//cout << " R2_side = " << R2_side_Correl << " +/- " << R2_side_Correl_err << endl;
		break;
	case 'l':
		R2_long[iKT][iKphi] = get_fit_results(1, solver_ptr)*hbarC*hbarC;
		R2_long_err[iKT][iKphi] = c*get_fit_err(1, covariance_ptr)*hbarC*hbarC;
		//cout << " R2_long = " << R2_long_Correl << " +/- " << R2_long_Correl_err << endl;
		break;
}

  //clean up
  gsl_matrix_free (covariance_ptr);
  gsl_rng_free (rng_ptr);

  delete[] Correlfun1D_data.q;
  delete[] Correlfun1D_data.y;
  delete[] Correlfun1D_data.sigma;

  gsl_multifit_fdfsolver_free (solver_ptr);  // free up the solver

  return;
}

//*********************************************************************
// 1D case
//*********************************************************************
//  Simple function to print results of each iteration in nice format
int CorrelationFunction::print_fit_state_1D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr)
{
  cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
  cout.precision (5);		// digits in doubles

  int width = 15;		// setw width for output
  cout << scientific
    << "iteration " << iteration << ": "
    << "x = {" << setw (width) << gsl_vector_get (solver_ptr->x, 0)
    << setw (width) << gsl_vector_get (solver_ptr->x, 1)
    << "}, |f(x)| = " << scientific << gsl_blas_dnrm2 (solver_ptr->f) 
    << endl << endl;

  return 0;
}
//*********************************************************************
//  Function returning the residuals for each point; that is, the 
//  difference of the fit function using the current parameters
//  and the data to be fit.
int Fittarget_correlfun1D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr)
{
  size_t n = ((struct Correlationfunction1D_data *) params_ptr)->data_length;
  double *q = ((struct Correlationfunction1D_data *) params_ptr)->q;
  double *y = ((struct Correlationfunction1D_data *) params_ptr)->y;
  double *sigma = ((struct Correlationfunction1D_data *) params_ptr)->sigma;

  //fit parameters
  double lambda = gsl_vector_get (xvec_ptr, 0);
  double R2 = gsl_vector_get (xvec_ptr, 1);

  size_t i;

  for (i = 0; i < n; i++)
  {
      double Yi = lambda*exp(- q[i]*q[i]*R2);
      gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
  }

  return GSL_SUCCESS;
}

//*********************************************************************
//  Function returning the Jacobian of the residual function
int Fittarget_correlfun1D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr)
{
  size_t n = ((struct Correlationfunction1D_data *) params_ptr)->data_length;
  double *q = ((struct Correlationfunction1D_data *) params_ptr)->q;
  double *sigma = ((struct Correlationfunction1D_data *) params_ptr)->sigma;

  //fit parameters
  double lambda = gsl_vector_get (xvec_ptr, 0);
  double R2 = gsl_vector_get (xvec_ptr, 1);

  size_t i;

  for (i = 0; i < n; i++)
  {
      // Jacobian matrix J(i,j) = dfi / dxj, 
      // where fi = (Yi - yi)/sigma[i],      
      //       Yi = A * exp(-lambda * i) + b 
      // and the xj are the parameters (A,lambda,b) 
      double sig = sigma[i];

      //derivatives
      double common_elemt = exp(- q[i]*q[i]*R2);
      
      gsl_matrix_set (Jacobian_ptr, i, 0, common_elemt/sig);
      //gsl_matrix_set (Jacobian_ptr, i, 1, - lambda*q[i]*q[i]*R2*common_elemt/sig);
      gsl_matrix_set (Jacobian_ptr, i, 1, - lambda*q[i]*q[i]*common_elemt/sig);
  }
  return GSL_SUCCESS;
}

//*********************************************************************
//  Function combining the residual function and its Jacobian
int Fittarget_correlfun1D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr)
{
  Fittarget_correlfun1D_f(xvec_ptr, params_ptr, f_ptr);
  Fittarget_correlfun1D_df(xvec_ptr, params_ptr, Jacobian_ptr);

  return GSL_SUCCESS;
}

//*********************************************************************
//  Function to return the i'th best-fit parameter
inline double CorrelationFunction::get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr)
{
  return gsl_vector_get (solver_ptr->x, i);
}

//*********************************************************************
//  Function to retrieve the square root of the diagonal elements of
//   the covariance matrix.
inline double CorrelationFunction::get_fit_err (int i, gsl_matrix * covariance_ptr)
{
  return sqrt (gsl_matrix_get (covariance_ptr, i, i));
}

//End of file
