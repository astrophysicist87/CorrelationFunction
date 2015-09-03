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
		Cal_correlationfunction_1D(iKT, iKphi);
		Fit_Correlationfunction1D('o', iKT, iKphi);
		Fit_Correlationfunction1D('s', iKT, iKphi);
		Fit_Correlationfunction1D('l', iKT, iKphi);
		Output_Correlationfunction_1D(iKT, iKphi, folderindex);
	 }
	 else if (corrfuncdim == 3)
	 {
		Cal_correlationfunction_3D(iKT, iKphi);
		if (lambdaflag)
			Fit_Correlationfunction3D_withlambda(iKT, iKphi);
		else Fit_Correlationfunction3D(iKT, iKphi);
		Output_Correlationfunction_3D(iKT, iKphi, folderindex);
	 }
      }
      R2_Fourier_transform(iKT, global_plane_psi);
   }
}

void CorrelationFunction::Cal_correlationfunction_1D(int iKT, int iKphi)
{
   if(fabs(K_y) > 1e-16)
   {
       cout<<"CorrelationFunction:: not support for y not equals 0 yet!" << endl;
       return;
   }
   
//int counter = 0;

   double mass = particle_mass;
   double local_K_T = K_T[iKT];
   double localK_phi = K_phi[iKphi];
   double cosK_phi = cos(localK_phi);
   double sinK_phi = sin(localK_phi);
   double error = 1e-4;

   //cout << "For Kt = " << 0. << " and Kphi = " << 0. << ", generating the 1d slices of the correlation function along q_out, q_side, and q_long direction... " << endl;
		  
   for(int i = 0; i < qnpts; i++)
   {
      //cout << "calculating q_mu = " << q_out[i] << endl;
      double values[3];
      for (int ops = 0; ops < 3; ops++)
         values[ops] = 0.0;
      for (int l = 0; l < 3; l++)
      //for (int l = 0; l < 2; l++)
      //for (int l = 0; l < 1; l++)
      {
         double local_q_out=0.0, local_q_side=0.0, local_q_long=0.0;
         switch (l)
         {
            case 0:
            {
               local_q_out  = q_out[i];
               local_q_side = 0.0e0;
               local_q_long = 0.0e0;
               break;
            }
   	      case 1:
            {
               local_q_out  = 0.0e0;
               local_q_side = q_side[i];
               local_q_long = 0.0e0;
               break;
            }
            case 2:
            {
               local_q_out  = 0.0e0;
               local_q_side = 0.0e0;
               local_q_long = q_long[i];
               break;
            }
            default:
            {
               cout << "error in assigning q values! "<< endl;
               break;
            }
         }

     	 double xsi  = local_K_T*local_K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
         double E1sq = xsi + local_K_T*local_q_out;
         double E2sq = xsi - local_K_T*local_q_out;
         double qt = sqrt(E1sq) - sqrt(E2sq);
         double qx = local_q_out*cosK_phi - local_q_side*sinK_phi;
         double qy = local_q_side*cosK_phi + local_q_out*sinK_phi;
         double qz = local_q_long;

         double integ1 = 0.0;                         //integ## are the numerators for the different q=const curves ("a" varies q_o, "b" varies q_s and "c" varies q_l) 
         double integ2 = 0.0;
spectra = 0.;

         for(int k = 0; k < Emissionfunction_length; k++)
         {
            double ss  = (*Emissionfunction_ptr)[k].data;
            double rpt = (*Emissionfunction_ptr)[k].r;
            double tpt = (*Emissionfunction_ptr)[k].t;
            double xpt = (*Emissionfunction_ptr)[k].x;
            double ypt = (*Emissionfunction_ptr)[k].y;
            double zpt = (*Emissionfunction_ptr)[k].z;

            for(int ii=0; ii<2; ii++)
            {
               zpt = zpt*(-1);   //using the symmetry along z axis
               double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
               integ1 += cos(arg)*ss;
               integ2 += sin(arg)*ss;
               spectra += ss;
            }
         }
         spectra /= (hbarC*hbarC*hbarC);
         integ1 = integ1/hbarC/hbarC/hbarC/spectra;
         integ2 = integ2/hbarC/hbarC/hbarC/spectra;
         double localvalue = integ1*integ1+integ2*integ2;
         values[l] = localvalue;
      }
      Correl_1D_out[i]  = values[0];
      Correl_1D_side[i] = values[1];
      Correl_1D_long[i] = values[2];
//cout << "CF: " << q_out[i] << "   " << Correl_1D_out[i] << endl;
      Correl_1D_out_err[i] = error;
      Correl_1D_side_err[i] = error;
      Correl_1D_long_err[i] = error;
   }
   return;
}

double CorrelationFunction::Emission_function(double p0, double pT, double pphi, double pz, FO_surf * FOcell_ptr, int local_pid)
{
	FO_surf FOcell = * FOcell_ptr;
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;
	
	double tau = FOcell.tau;
	double r = FOcell.r;
	double sin_phi = FOcell.sin_phi;
	double cos_phi = FOcell.cos_phi;

	double da0 = FOcell.da0;
	double da1 = FOcell.da1;
	double da2 = FOcell.da2;

	double vx = FOcell.vx;
	double vy = FOcell.vy;
	double gammaT = FOcell.gammaT;

	double pi00 = FOcell.pi00;
	double pi01 = FOcell.pi01;
	double pi02 = FOcell.pi02;
	double pi11 = FOcell.pi11;
	double pi12 = FOcell.pi12;
	double pi22 = FOcell.pi22;
	double pi33 = FOcell.pi33;

	double Tdec = FOcell.Tdec;
	double Pdec = FOcell.Pdec;
	double Edec = FOcell.Edec;

	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double deltaf_prefactor = 0.;
	if (use_delta_f) deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	pT = SPinterp_pT[ipt];
	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;
	
	//now get distribution function, emission function, etc.
	double f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
	//viscous corrections
	double deltaf = 0.;
	if (use_delta_f)
		deltaf = deltaf_prefactor * (1. - sign*f0)
			* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);
	//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
	double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);
	//ignore points where delta f is large or emission function goes negative from pdsigma
	if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
		S_p = 0.0;

	//S_p_withweight = S_p*tau*eta_s_weight[ieta];

	return (S_p);
}

void CorrelationFunction::Cal_correlationfunction_3D(int iKT, int iKphi)
{
   if(fabs(K_y) > 1e-16)
   {
       //cout<<"not support for y not equals 0 yet!" << endl;
       return;
   }

   double mass = particle_mass;
   double local_K_T = K_T[iKT];
   double localK_phi = K_phi[iKphi];
   double cosK_phi = cos(localK_phi);
   double sinK_phi = sin(localK_phi);
   double error = 1e-4;

   //cout << "generating correlation function in 3D... " << endl;
   for(int i = 0; i < qnpts; i++)  // q_out loop
   {
      double local_q_out = q_out[i];
      //cout << "q_out = " << local_q_out << endl;
      for(int j = 0; j < qnpts; j++)  // q_side loop
      {
         double local_q_side = q_side[j];
         for(int k = 0; k < qnpts; k++)  // q_long loop
         {
            double local_q_long = q_long[k];
            double integ1 = 0.0;                         
            double integ2 = 0.0;
            double sum = 0.0;

     	    double xsi = local_K_T*local_K_T + mass*mass + (local_q_out*local_q_out + local_q_side*local_q_side + local_q_long*local_q_long)/4.0;  //Set Xsi
            double E1sq = xsi + local_K_T*local_q_out;
            double E2sq = xsi - local_K_T*local_q_out;
            double qt = sqrt(E1sq) - sqrt(E2sq);
            double qx = local_q_out*cosK_phi - local_q_side*sinK_phi;
            double qy = local_q_side*cosK_phi + local_q_out*sinK_phi;
            double qz = local_q_long;
            
            for(int m = 0; m < Emissionfunction_length; m++)
            {
            double ss  = (*Emissionfunction_ptr)[k].data;
            double tpt = (*Emissionfunction_ptr)[k].t;
            double xpt = (*Emissionfunction_ptr)[k].x;
            double ypt = (*Emissionfunction_ptr)[k].y;
            double zpt = (*Emissionfunction_ptr)[k].z;

               for(int ii=0; ii<2; ii++)
               {
                  zpt = zpt*(-1);
               
                  double arg = (tpt*qt - (qx*xpt + qy*ypt + qz*zpt))/hbarC;
                  integ1 += cos(arg)*ss;
                  integ2 += sin(arg)*ss;
               }
            }
            integ1 = integ1/hbarC/hbarC/hbarC/spectra;
            integ2 = integ2/hbarC/hbarC/hbarC/spectra;
            sum = integ1*integ1+integ2*integ2;
            Correl_3D[i][j][k] = sum;
            Correl_3D_err[i][j][k] = error;
         }
      }
   }
   return;
}

int CorrelationFunction::Read_correlationfunction_3D(int iKT, int iKphi)
{
   if(fabs(K_y) > 1e-16)
   {
       //cout<<"not support for y not equals 0 yet!" << endl;
       return (0);
   }

  ostringstream corrfn_stream;
  corrfn_stream << "output/stored_results/results_exact_kt_" << K_T << "/correlfunct3D" << "_" << particle_name << "_kt_" << K_T << "_phi_" << K_phi << ".dat";
if ( !fexists( corrfn_stream.str().c_str() ) )
{
	cerr << corrfn_stream.str().c_str() << ": File not found" << endl;
	return (1);
}
  ifstream corrfn(corrfn_stream.str().c_str());

   for(int i = 0; i < qnpts; i++)  // q_out loop
   {
      for(int j = 0; j < qnpts; j++)  // q_side loop
      {
         for(int k = 0; k < qnpts; k++)  // q_long loop
         {
            corrfn >> q_out[i];
            corrfn >> q_side[j];
            corrfn >> q_long[k];
            corrfn >> Correl_3D[i][j][k];
            corrfn >> Correl_3D_err[i][j][k];
         }
      }
   }

   return (0);
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
//corrfn_stream << path;
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

void CorrelationFunction::Fit_Correlationfunction3D(int iKT, int iKphi)
{
  const size_t data_length = qnpts*qnpts*qnpts;  // # of points
  const size_t n_para = 4;  // # of parameters

  // allocate space for a covariance matrix of size p by p
  gsl_matrix *covariance_ptr = gsl_matrix_alloc (n_para, n_para);

  // allocate and setup for generating gaussian distibuted random numbers
  gsl_rng_env_setup ();
  const gsl_rng_type *type = gsl_rng_default;
  gsl_rng *rng_ptr = gsl_rng_alloc (type);

  //set up test data
  struct Correlationfunction3D_data Correlfun3D_data;
  Correlfun3D_data.data_length = data_length;
  Correlfun3D_data.q_o = new double [data_length];
  Correlfun3D_data.q_s = new double [data_length];
  Correlfun3D_data.q_l = new double [data_length];
  Correlfun3D_data.y = new double [data_length];
  Correlfun3D_data.sigma = new double [data_length];

  int idx = 0;
  for(int i=0; i<qnpts; i++)
  {
    for(int j=0; j<qnpts; j++)
    {
      for(int k=0; k<qnpts; k++)
      {
         Correlfun3D_data.q_o[idx] = q_out[i];
         Correlfun3D_data.q_s[idx] = q_side[j];
         Correlfun3D_data.q_l[idx] = q_long[k];
         // This sets up the data to be fitted, with gaussian noise added
         // Correlfun3D_data.y[idx] = 1.0*exp( - 0.81*q_out[i]*q_out[i] - 1.21*q_side[j]*q_side[j] - 4.0*q_long[k]*q_long[k] - 0.25*q_out[i]*q_side[j]) + gsl_ran_gaussian(rng_ptr, error);
         Correlfun3D_data.y[idx] = Correl_3D[i][j][k];
         Correlfun3D_data.sigma[idx] = Correl_3D_err[i][j][k];
         //Correlfun3D_data.sigma[idx] = 1e-2;
         idx++;
      }
    }
  }

  double para_init[n_para] = { 1.0, 1.0, 1.0, 1.0 };  // initial guesse of parameters

  gsl_vector_view xvec_ptr = gsl_vector_view_array (para_init, n_para);
  
  // set up the function to be fit 
  gsl_multifit_function_fdf target_func;
  target_func.f = &Fittarget_correlfun3D_f;        // the function of residuals
  target_func.df = &Fittarget_correlfun3D_df;      // the gradient of this function
  target_func.fdf = &Fittarget_correlfun3D_fdf;    // combined function and gradient
  target_func.n = data_length;              // number of points in the data set
  target_func.p = n_para;              // number of parameters in the fit function
  target_func.params = &Correlfun3D_data;  // structure with the data and error bars

  const gsl_multifit_fdfsolver_type *type_ptr = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver *solver_ptr 
       = gsl_multifit_fdfsolver_alloc (type_ptr, data_length, n_para);
  gsl_multifit_fdfsolver_set (solver_ptr, &target_func, &xvec_ptr.vector);

  size_t iteration = 0;         // initialize iteration counter
  //print_fit_state_3D (iteration, solver_ptr);
  int status;  		// return value from gsl function calls (e.g., error)
  do
  {
      iteration++;
      
      // perform a single iteration of the fitting routine
      status = gsl_multifit_fdfsolver_iterate (solver_ptr);

      // print out the status of the fit
      //cout << "status = " << gsl_strerror (status) << endl;

      // customized routine to print out current parameters
      //print_fit_state_3D (iteration, solver_ptr);

      if (status)    // check for a nonzero status code
      {
          break;  // this should only happen if an error code is returned 
      }

      // test for convergence with an absolute and relative error (see manual)
      status = gsl_multifit_test_delta (solver_ptr->dx, solver_ptr->x, 
                                        fit_tolarence, fit_tolarence);
  }
  while (status == GSL_CONTINUE && iteration < fit_max_iterations);

//cerr >> "iteration = " << iteration << endl;

  // calculate the covariance matrix of the best-fit parameters
  gsl_multifit_covar (solver_ptr->J, 0.0, covariance_ptr);

  // print out the covariance matrix using the gsl function (not elegant!)
  //cout << endl << "Covariance matrix: " << endl;
  //gsl_matrix_fprintf (stdout, covariance_ptr, "%g");

  //cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
  //cout.precision (5);		                // # of digits in doubles

  int width = 7;		// setw width for output
  //cout << endl << "Best fit results:" << endl;
  //cout << "R2o = " << setw (width) << get_fit_results (0, solver_ptr)
  //  << " +/- " << setw (width) << get_fit_err (0, covariance_ptr) << endl;

  //cout << "R2s      = " << setw (width) << get_fit_results (1, solver_ptr)
  //  << " +/- " << setw (width) << get_fit_err (1, covariance_ptr) << endl;

  //cout << "R2l      = " << setw (width) << get_fit_results (2, solver_ptr)
  //  << " +/- " << setw (width) << get_fit_err (2, covariance_ptr) << endl;
  
  //cout << "R2os      = " << setw (width) << get_fit_results (3, solver_ptr)
  //  << " +/- " << setw (width) << get_fit_err (3, covariance_ptr) << endl;
    
  //cout << "status = " << gsl_strerror (status) << endl;
  //cout << "--------------------------------------------------------------------" << endl;

  double chi = gsl_blas_dnrm2(solver_ptr->f);
  double dof = data_length - n_para;
  double c = GSL_MAX_DBL(1, chi/sqrt(dof));

  lambda_Correl[iKT][iKphi] = 1.0;
  lambda_Correl_err[iKT][iKphi] = 0.0;
  R2_out[iKT][iKphi] = fabs(get_fit_results(0, solver_ptr))*hbarC*hbarC;
  R2_side[iKT][iKphi] = fabs(get_fit_results(1, solver_ptr))*hbarC*hbarC;
  R2_long[iKT][iKphi] = fabs(get_fit_results(2, solver_ptr))*hbarC*hbarC;
  R2_outside[iKT][iKphi] = fabs(get_fit_results(3, solver_ptr))*hbarC*hbarC;
  R2_out_err[iKT][iKphi] = c*get_fit_err(0, covariance_ptr)*hbarC*hbarC;
  R2_side_err[iKT][iKphi] = c*get_fit_err(1, covariance_ptr)*hbarC*hbarC;
  R2_long_err[iKT][iKphi] = c*get_fit_err(2, covariance_ptr)*hbarC*hbarC;
  R2_outside_err[iKT][iKphi] = c*get_fit_err(3, covariance_ptr)*hbarC*hbarC;

  //cout << "final results: " << endl;
  //cout << scientific << setw(10) << setprecision(5) 
  //     << "chisq/dof = " << chi*chi/dof << endl;
  //cout << scientific << setw(10) << setprecision(5) 
  //     << " lambda = " << lambda_CorrelGF << " +/- " << lambda_CorrelGF_err << endl;
  //cout << " R2_out = " << R2_out_CorrelGF << " +/- " << R2_out_CorrelGF_err << endl;
  //cout << " R2_side = " << R2_side_CorrelGF << " +/- " << R2_side_CorrelGF_err << endl;
  //cout << " R2_long = " << R2_long_CorrelGF << " +/- " << R2_long_CorrelGF_err << endl;
  //cout << " R2_outside = " << R2_outside_CorrelGF << " +/- " << R2_outside_CorrelGF_err << endl;

  //clean up
  gsl_matrix_free (covariance_ptr);
  gsl_rng_free (rng_ptr);

  delete[] Correlfun3D_data.q_o;
  delete[] Correlfun3D_data.q_s;
  delete[] Correlfun3D_data.q_l;
  delete[] Correlfun3D_data.y;
  delete[] Correlfun3D_data.sigma;

  gsl_multifit_fdfsolver_free (solver_ptr);  // free up the solver

  return;
}

void CorrelationFunction::Fit_Correlationfunction3D_withlambda(int iKT, int iKphi)
{
  const size_t data_length = qnpts*qnpts*qnpts;  // # of points
  const size_t n_para = 5;  // # of parameters

  // allocate space for a covariance matrix of size p by p
  gsl_matrix *covariance_ptr = gsl_matrix_alloc (n_para, n_para);

  // allocate and setup for generating gaussian distibuted random numbers
  gsl_rng_env_setup ();
  const gsl_rng_type *type = gsl_rng_default;
  gsl_rng *rng_ptr = gsl_rng_alloc (type);

  //set up data
  struct Correlationfunction3D_data Correlfun3D_data;
  Correlfun3D_data.data_length = data_length;
  Correlfun3D_data.q_o = new double [data_length];
  Correlfun3D_data.q_s = new double [data_length];
  Correlfun3D_data.q_l = new double [data_length];
  Correlfun3D_data.y = new double [data_length];
  Correlfun3D_data.sigma = new double [data_length];

  int idx = 0;
  for(int i=0; i<qnpts; i++)
  {
    for(int j=0; j<qnpts; j++)
    {
      for(int k=0; k<qnpts; k++)
      {
         Correlfun3D_data.q_o[idx] = q_out[i];
         Correlfun3D_data.q_s[idx] = q_side[j];
         Correlfun3D_data.q_l[idx] = q_long[k];
         Correlfun3D_data.y[idx] = Correl_3D[i][j][k];
         Correlfun3D_data.sigma[idx] = Correl_3D_err[i][j][k];
         idx++;
      }
    }
  }

  double para_init[n_para] = { 1.0, 1.0, 1.0, 1.0, 1.0 };  // initial guesse of parameters

  gsl_vector_view xvec_ptr = gsl_vector_view_array (para_init, n_para);
  
  // set up the function to be fit 
  gsl_multifit_function_fdf target_func;
  target_func.f = &Fittarget_correlfun3D_f_withlambda;        // the function of residuals
  target_func.df = &Fittarget_correlfun3D_df_withlambda;      // the gradient of this function
  target_func.fdf = &Fittarget_correlfun3D_fdf_withlambda;    // combined function and gradient
  target_func.n = data_length;              // number of points in the data set
  target_func.p = n_para;              // number of parameters in the fit function
  target_func.params = &Correlfun3D_data;  // structure with the data and error bars

  const gsl_multifit_fdfsolver_type *type_ptr = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver *solver_ptr 
       = gsl_multifit_fdfsolver_alloc (type_ptr, data_length, n_para);
  gsl_multifit_fdfsolver_set (solver_ptr, &target_func, &xvec_ptr.vector);

  size_t iteration = 0;         // initialize iteration counter
  //print_fit_state_3D (iteration, solver_ptr);
  int status;  		// return value from gsl function calls (e.g., error)
  do
  {
      iteration++;
      
      // perform a single iteration of the fitting routine
      status = gsl_multifit_fdfsolver_iterate (solver_ptr);

      // print out the status of the fit
      //cout << "status = " << gsl_strerror (status) << endl;

      // customized routine to print out current parameters
      //print_fit_state_3D (iteration, solver_ptr);

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
  //cout << "lambda      = " << setw (width) << get_fit_results (0, solver_ptr)
  //  << " +/- " << setw (width) << get_fit_err (0, covariance_ptr) << endl;

  //cout << "R2o = " << setw (width) << get_fit_results (1, solver_ptr)
  //  << " +/- " << setw (width) << get_fit_err (1, covariance_ptr) << endl;

  //cout << "R2s      = " << setw (width) << get_fit_results (2, solver_ptr)
  //  << " +/- " << setw (width) << get_fit_err (2, covariance_ptr) << endl;

  //cout << "R2l      = " << setw (width) << get_fit_results (3, solver_ptr)
  //  << " +/- " << setw (width) << get_fit_err (3, covariance_ptr) << endl;
  
  //cout << "R2os      = " << setw (width) << get_fit_results (4, solver_ptr)
  //  << " +/- " << setw (width) << get_fit_err (4, covariance_ptr) << endl;
    
  //cout << "status = " << gsl_strerror (status) << endl;
  //cout << "--------------------------------------------------------------------" << endl;

  double chi = gsl_blas_dnrm2(solver_ptr->f);
  double dof = data_length - n_para;
  double c = GSL_MAX_DBL(1, chi/sqrt(dof));

  lambda_Correl[iKT][iKphi] = get_fit_results(0, solver_ptr);
  R2_out[iKT][iKphi] = fabs(get_fit_results(1, solver_ptr))*hbarC*hbarC;
  R2_side[iKT][iKphi] = fabs(get_fit_results(2, solver_ptr))*hbarC*hbarC;
  R2_long[iKT][iKphi] = fabs(get_fit_results(3, solver_ptr))*hbarC*hbarC;
  R2_outside[iKT][iKphi] = fabs(get_fit_results(4, solver_ptr))*hbarC*hbarC;
  lambda_Correl_err[iKT][iKphi] = c*get_fit_err(0, covariance_ptr);
  R2_out_err[iKT][iKphi] = c*get_fit_err(1, covariance_ptr)*hbarC*hbarC;
  R2_side_err[iKT][iKphi] = c*get_fit_err(2, covariance_ptr)*hbarC*hbarC;
  R2_long_err[iKT][iKphi] = c*get_fit_err(3, covariance_ptr)*hbarC*hbarC;
  R2_outside_err[iKT][iKphi] = c*get_fit_err(4, covariance_ptr)*hbarC*hbarC;

  //cout << "final results: " << endl;
  //cout << scientific << setw(10) << setprecision(5) 
  //     << "chisq/dof = " << chi*chi/dof << endl;
  //cout << scientific << setw(10) << setprecision(5) 
  //     << " lambda = " << lambda_CorrelGF << " +/- " << lambda_CorrelGF_err << endl;
  //cout << " R2_out = " << R2_out_CorrelGF << " +/- " << R2_out_CorrelGF_err << endl;
  //cout << " R2_side = " << R2_side_CorrelGF << " +/- " << R2_side_CorrelGF_err << endl;
  //cout << " R2_long = " << R2_long_CorrelGF << " +/- " << R2_long_CorrelGF_err << endl;
  //cout << " R2_outside = " << R2_outside_CorrelGF << " +/- " << R2_outside_CorrelGF_err << endl;

  //clean up
  gsl_matrix_free (covariance_ptr);
  gsl_rng_free (rng_ptr);

  delete[] Correlfun3D_data.q_o;
  delete[] Correlfun3D_data.q_s;
  delete[] Correlfun3D_data.q_l;
  delete[] Correlfun3D_data.y;
  delete[] Correlfun3D_data.sigma;

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
// 3D case
//*********************************************************************
//  Simple function to print results of each iteration in nice format
int CorrelationFunction::print_fit_state_3D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr)
{
  cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
  cout.precision (5);		// digits in doubles

  int width = 15;		// setw width for output
  cout << scientific
    << "iteration " << iteration << ": "
    << "  x = {" << setw (width) << gsl_vector_get (solver_ptr->x, 0)
    << setw (width) << gsl_vector_get (solver_ptr->x, 1)
    << setw (width) << gsl_vector_get (solver_ptr->x, 2)
    << setw (width) << gsl_vector_get (solver_ptr->x, 3)
    << "}, |f(x)| = " << scientific << gsl_blas_dnrm2 (solver_ptr->f) 
    << endl << endl;

  return 0;
}
//  Simple function to print results of each iteration in nice format
int CorrelationFunction::print_fit_state_3D_withlambda (size_t iteration, gsl_multifit_fdfsolver * solver_ptr)
{
  cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
  cout.precision (5);		// digits in doubles

  int width = 15;		// setw width for output
  cout << scientific
    << "iteration " << iteration << ": "
    << "  x = {" << setw (width) << gsl_vector_get (solver_ptr->x, 0)
    << setw (width) << gsl_vector_get (solver_ptr->x, 1)
    << setw (width) << gsl_vector_get (solver_ptr->x, 2)
    << setw (width) << gsl_vector_get (solver_ptr->x, 3)
    << setw (width) << gsl_vector_get (solver_ptr->x, 4)
    << "}, |f(x)| = " << scientific << gsl_blas_dnrm2 (solver_ptr->f) 
    << endl << endl;

  return 0;
}
//*********************************************************************
//  Function returning the residuals for each point; that is, the 
//  difference of the fit function using the current parameters
//  and the data to be fit.
int Fittarget_correlfun3D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr)
{
  size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
  double *q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
  double *q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
  double *q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
  double *y = ((struct Correlationfunction3D_data *) params_ptr)->y;
  double *sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

  //fit parameters
  double R2_o = gsl_vector_get (xvec_ptr, 0);
  double R2_s = gsl_vector_get (xvec_ptr, 1);
  double R2_l = gsl_vector_get (xvec_ptr, 2);
  double R2_os = gsl_vector_get (xvec_ptr, 3);

  size_t i;

  for (i = 0; i < n; i++)
  {
      double Yi = exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s
                   - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
      gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
  }

  return GSL_SUCCESS;
}

int Fittarget_correlfun3D_f_withlambda (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr)
{
  size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
  double *q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
  double *q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
  double *q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
  double *y = ((struct Correlationfunction3D_data *) params_ptr)->y;
  double *sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

  //fit parameters
  double lambda = gsl_vector_get (xvec_ptr, 0);
  double R2_o = gsl_vector_get (xvec_ptr, 1);
  double R2_s = gsl_vector_get (xvec_ptr, 2);
  double R2_l = gsl_vector_get (xvec_ptr, 3);
  double R2_os = gsl_vector_get (xvec_ptr, 4);

  size_t i;

  for (i = 0; i < n; i++)
  {
      //double Yi = lambda*exp(- q_l[i]*q_l[i]*R_l*R_l - q_s[i]*q_s[i]*R_s*R_s
      //             - q_o[i]*q_o[i]*R_o*R_o - q_o[i]*q_s[i]*R_os*R_os);
      double Yi = lambda*exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s
                   - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
      gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
  }

  return GSL_SUCCESS;
}

//*********************************************************************
//  Function returning the Jacobian of the residual function
int Fittarget_correlfun3D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr)
{
  size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
  double *q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
  double *q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
  double *q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
  double *sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

  //fit parameters
  double R2_o = gsl_vector_get (xvec_ptr, 0);
  double R2_s = gsl_vector_get (xvec_ptr, 1);
  double R2_l = gsl_vector_get (xvec_ptr, 2);
  double R2_os = gsl_vector_get (xvec_ptr, 3);

  size_t i;

  for (i = 0; i < n; i++)
  {
      // Jacobian matrix J(i,j) = dfi / dxj, 
      // where fi = (Yi - yi)/sigma[i],      
      //       Yi = A * exp(-lambda * i) + b 
      // and the xj are the parameters (A,lambda,b) 
      double sig = sigma[i];

      //derivatives
      //double common_elemt = exp(- q_l[i]*q_l[i]*R_l*R_l - q_s[i]*q_s[i]*R_s*R_s
      //             - q_o[i]*q_o[i]*R_o*R_o - q_o[i]*q_s[i]*R_os*R_os);
      double common_elemt = exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s
                   - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
      
      gsl_matrix_set (Jacobian_ptr, i, 0, - q_o[i]*q_o[i]*common_elemt/sig);
      gsl_matrix_set (Jacobian_ptr, i, 1, - q_s[i]*q_s[i]*common_elemt/sig);
      gsl_matrix_set (Jacobian_ptr, i, 2, - q_l[i]*q_l[i]*common_elemt/sig);
      gsl_matrix_set (Jacobian_ptr, i, 3, - 2.*q_o[i]*q_s[i]*common_elemt/sig);
  }
  return GSL_SUCCESS;
}

int Fittarget_correlfun3D_df_withlambda (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr)
{
  size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
  double *q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
  double *q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
  double *q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
  double *sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

  //fit parameters
  double lambda = gsl_vector_get (xvec_ptr, 0);
  double R2_o = gsl_vector_get (xvec_ptr, 1);
  double R2_s = gsl_vector_get (xvec_ptr, 2);
  double R2_l = gsl_vector_get (xvec_ptr, 3);
  double R2_os = gsl_vector_get (xvec_ptr, 4);

  size_t i;

  for (i = 0; i < n; i++)
  {
      // Jacobian matrix J(i,j) = dfi / dxj, 
      // where fi = (Yi - yi)/sigma[i],      
      //       Yi = A * exp(-lambda * i) + b 
      // and the xj are the parameters (A,lambda,b) 
      double sig = sigma[i];

      //derivatives
      double common_elemt = exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s
                   - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
      
      gsl_matrix_set (Jacobian_ptr, i, 0, common_elemt/sig);
      gsl_matrix_set (Jacobian_ptr, i, 1, - lambda*q_o[i]*q_o[i]*common_elemt/sig);
      gsl_matrix_set (Jacobian_ptr, i, 2, - lambda*q_s[i]*q_s[i]*common_elemt/sig);
      gsl_matrix_set (Jacobian_ptr, i, 3, - lambda*q_l[i]*q_l[i]*common_elemt/sig);
      gsl_matrix_set (Jacobian_ptr, i, 4, - 2.*lambda*q_o[i]*q_s[i]*common_elemt/sig);
  }
  return GSL_SUCCESS;
}

//*********************************************************************
//  Function combining the residual function and its Jacobian
int Fittarget_correlfun3D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr)
{
  Fittarget_correlfun3D_f(xvec_ptr, params_ptr, f_ptr);
  Fittarget_correlfun3D_df(xvec_ptr, params_ptr, Jacobian_ptr);

  return GSL_SUCCESS;
}
int Fittarget_correlfun3D_fdf_withlambda (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr)
{
  Fittarget_correlfun3D_f_withlambda(xvec_ptr, params_ptr, f_ptr);
  Fittarget_correlfun3D_df_withlambda(xvec_ptr, params_ptr, Jacobian_ptr);

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
