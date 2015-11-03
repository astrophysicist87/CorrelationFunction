#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>

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

double CorrelationFunction::get_Q()
{
	double smin = (m2+m3)*(m2+m3);
	double smax = (Mres-mass)*(Mres-mass);
	double sum = 0.;
	
	for (int is = 0; is < n_s_pts; ++is)
	{
		double sp = NEW_s_pts[is];
		double f1 = (Mres+mass)*(Mres+mass) - sp;
		double f2 = smax - sp;
		double f3 = smin - sp;
		double f4 = (m2-m3)*(m2-m3) - sp;
		sum += NEW_s_wts[is]*sqrt(f1*f2*f3*f4)/(sp+1.e-15);
	}

	return sum;
}

double CorrelationFunction::g(double s)
{
	double gs_pstar_loc = sqrt( ((Mres+mass)*(Mres+mass) - s)*((Mres-mass)*(Mres-mass) - s) )/(2.0*Mres);
	double g_res = br/(4.*M_PI*gs_pstar_loc);
	if (n_body == 3 || n_body == 4)		//both set up to work the same way
	{
		double pre_f = (Mres*br)/(2.*M_PI*s);
		double num = sqrt( (s - (m2+m3)*(m2+m3)) * (s - (m2-m3)*(m2-m3)) );
		double den = Qfunc;
		g_res = pre_f * num / den;
	}

	return g_res;
}

void CorrelationFunction::Do_resonance_integrals(int parent_resonance_particle_id, int daughter_particle_id, int decay_channel)
{
	time_t rawtime;
  	struct tm * timeinfo;

	Allocate_resonance_running_sum_vectors();

	res_sign_info = sign_of_dN_dypTdpTdphi_moments[parent_resonance_particle_id];
	res_log_info = ln_dN_dypTdpTdphi_moments[parent_resonance_particle_id];
	res_moments_info = dN_dypTdpTdphi_moments[parent_resonance_particle_id];

	n_body = current_reso_nbody;

	if (n_body == 2)
	{
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			double local_pT = SPinterp_pT[ipt];
			double local_pphi = SPinterp_pphi[ipphi];
			Zero_resonance_running_sum_vector(ssum_vec);
			Zero_resonance_running_sum_vector(vsum_vec);
			Zero_resonance_running_sum_vector(zetasum_vec);
			Zero_resonance_running_sum_vector(Csum_vec);
			Load_decay_channel_info(decay_channel, local_pT, local_pphi);	// set decay channel information

			//then g(s) is delta-function, skip s-integration entirely
			//double s_loc = m2*m2;
			for (int iv = 0; iv < n_v_pts; ++iv)
			{
				time (&rawtime);
				timeinfo = localtime (&rawtime);
				Zero_resonance_running_sum_vector(zetasum_vec);
				for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
				{
					Zero_resonance_running_sum_vector(Csum_vec);
					double PKT = VEC_n2_PT[iv][izeta];
					double PKY = VEC_n2_P_Y[iv];
					double PKphi = VEC_n2_PPhi_tilde[iv][izeta];
					for (int tempidx = 1; tempidx <= 2; ++tempidx)
					{
						if (tempidx != 1)
							PKphi = VEC_n2_PPhi_tildeFLIP[iv][izeta];		//also takes Pp --> Pm
						Edndp3(PKT, PKphi, Csum_vec);
					}												// end of tempidx sum
					for (int iq = 0; iq < qnpts; ++iq)
					for (int iqax = 0; iqax < 3; ++iqax)
					for (int itrig = 0; itrig < 2; ++itrig)
						zetasum_vec[iq][iqax][itrig] += VEC_n2_zeta_factor[iv][izeta]*Csum_vec[iq][iqax][itrig];
				}													// end of zeta sum
				for (int iq = 0; iq < qnpts; ++iq)
				for (int iqax = 0; iqax < 3; ++iqax)
				for (int itrig = 0; itrig < 2; ++itrig)
					vsum_vec[iq][iqax][itrig] += VEC_n2_v_factor[iv]*zetasum_vec[iq][iqax][itrig];
			}														// end of v sum
			for (int iq = 0; iq < qnpts; ++iq)
			for (int iqax = 0; iqax < 3; ++iqax)
			for (int itrig = 0; itrig < 2; ++itrig)
				ssum_vec[iq][iqax][itrig] += Mres*VEC_n2_s_factor*vsum_vec[iq][iqax][itrig];

			//update all gridpoints for all daughter moments
			for (int iq = 0; iq < qnpts; ++iq)
			for (int iqax = 0; iqax < 3; ++iqax)
			for (int itrig = 0; itrig < 2; ++itrig)
			{
				dN_dypTdpTdphi_moments[daughter_particle_id][iq][iqax][itrig][ipt][ipphi] += ssum_vec[iq][iqax][itrig] / fraction_of_resonances;
				double temp = dN_dypTdpTdphi_moments[daughter_particle_id][iq][iqax][itrig][ipt][ipphi];
				ln_dN_dypTdpTdphi_moments[daughter_particle_id][iq][iqax][itrig][ipt][ipphi] = log(abs(temp));
				sign_of_dN_dypTdpTdphi_moments[daughter_particle_id][iq][iqax][itrig][ipt][ipphi] = sgn(temp);
			}
	
			if (isnan(dN_dypTdpTdphi_moments[daughter_particle_id][0][0][0][ipt][ipphi]
					+ dN_dypTdpTdphi_moments[daughter_particle_id][0][0][1][ipt][ipphi]))
			{
				*global_out_stream_ptr << "ERROR: NaNs encountered!" << endl
										<< "dN_dypTdpTdphi_moments[daughter_particle_id][0][0][0][" << ipt << "][" << ipphi << "] = "
										<< setw(8) << setprecision(15)
										<< dN_dypTdpTdphi_moments[daughter_particle_id][0][0][0][ipt][ipphi] << endl
										<< "dN_dypTdpTdphi_moments[daughter_particle_id][0][0][1][" << ipt << "][" << ipphi << "] = "
										<< setw(8) << setprecision(15)
										<< dN_dypTdpTdphi_moments[daughter_particle_id][0][0][1][ipt][ipphi] << endl
										<< "  --> pt = " << local_pT << std::endl
										<< "  --> pphi = " << local_pphi << std::endl
										<< "daughter_particle_id = " << daughter_particle_id << endl
										<< "parent_resonance_particle_id = " << parent_resonance_particle_id << endl
										<< "  --> Qfunc = " << Qfunc << endl
										<< "  --> n_body = " << n_body << endl
										<< "  --> gRES = " << gRES << endl
										<< "  --> Mres = " << Mres << endl
										<< "  --> mass = " << mass << endl
										<< "  --> Gamma = " << Gamma << endl
										<< "  --> br = " << br << endl
										<< "  --> m2 = " << m2 << endl
										<< "  --> m3 = " << m3 << endl << endl;
				exit(1);
			}
		}											// end of pT, pphi loops
	}												// end of nbody == 2
	else
	{
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			double local_pT = SPinterp_pT[ipt];
			double local_pphi = SPinterp_pphi[ipphi];
			Zero_resonance_running_sum_vector(ssum_vec);
			Zero_resonance_running_sum_vector(vsum_vec);
			Zero_resonance_running_sum_vector(zetasum_vec);
			Zero_resonance_running_sum_vector(Csum_vec);
			Load_decay_channel_info(decay_channel, local_pT, local_pphi);	// set decay channel information

			for (int is = 0; is < n_s_pts; ++is)
			{
				double vsum = 0.0;
 		  		Zero_resonance_running_sum_vector(vsum_vec);
				for (int iv = 0; iv < n_v_pts; ++iv)
				{
					Zero_resonance_running_sum_vector(zetasum_vec);
					for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
					{
						Zero_resonance_running_sum_vector(Csum_vec);
						double PKT = VEC_PT[is][iv][izeta];
						double PKY = VEC_P_Y[is][iv];
						double PKphi = VEC_PPhi_tilde[is][iv][izeta];
						for (int tempidx = 1; tempidx <= 2; ++tempidx)
						{
							if (tempidx != 1)
								PKphi = VEC_PPhi_tildeFLIP[is][iv][izeta];		//also takes Pp --> Pm
							Edndp3(PKT, PKphi, Csum_vec);
						}										// end of tempidx sum
                        for (int iq = 0; iq < qnpts; ++iq)
                        for (int iqax = 0; iqax < 3; ++iqax)
                        for (int itrig = 0; itrig < 2; ++itrig)
							zetasum_vec[iq][iqax][itrig] += VEC_zeta_factor[is][iv][izeta]*Csum_vec[iq][iqax][itrig];
					}											// end of zeta sum
                    for (int iq = 0; iq < qnpts; ++iq)
                    for (int iqax = 0; iqax < 3; ++iqax)
                    for (int itrig = 0; itrig < 2; ++itrig)
					    vsum_vec[iq][iqax][itrig] += VEC_v_factor[is][iv]*zetasum_vec[iq][iqax][itrig];
				}												// end of v sum
				for (int iq = 0; iq < qnpts; ++iq)
				for (int iqax = 0; iqax < 3; ++iqax)
				for (int itrig = 0; itrig < 2; ++itrig)
					ssum_vec[iq][iqax][itrig] += Mres*VEC_s_factor[is]*vsum_vec[iq][iqax][itrig];
			}													// end of s sum
			//update all gridpoints for daughter moments
			for (int iq = 0; iq < qnpts; ++iq)
			for (int iqax = 0; iqax < 3; ++iqax)
			for (int itrig = 0; itrig < 2; ++itrig)
			{
				dN_dypTdpTdphi_moments[daughter_particle_id][iq][iqax][itrig][ipt][ipphi] += ssum_vec[iq][iqax][itrig] / fraction_of_resonances;
				double temp = dN_dypTdpTdphi_moments[daughter_particle_id][iq][iqax][itrig][ipt][ipphi];
				ln_dN_dypTdpTdphi_moments[daughter_particle_id][iq][iqax][itrig][ipt][ipphi] = log(abs(temp));
				sign_of_dN_dypTdpTdphi_moments[daughter_particle_id][iq][iqax][itrig][ipt][ipphi] = sgn(temp);
			}
	
			if (isnan(dN_dypTdpTdphi_moments[daughter_particle_id][0][0][0][ipt][ipphi]
					+ dN_dypTdpTdphi_moments[daughter_particle_id][0][0][1][ipt][ipphi]))
			{
				*global_out_stream_ptr << "ERROR: NaNs encountered!" << endl
										<< "dN_dypTdpTdphi_moments[daughter_particle_id][0][0][0][" << ipt << "][" << ipphi << "] = "
										<< setw(8) << setprecision(15)
										<< dN_dypTdpTdphi_moments[daughter_particle_id][0][0][0][ipt][ipphi] << endl
										<< "dN_dypTdpTdphi_moments[daughter_particle_id][0][0][1][" << ipt << "][" << ipphi << "] = "
										<< setw(8) << setprecision(15)
										<< dN_dypTdpTdphi_moments[daughter_particle_id][0][0][1][ipt][ipphi] << endl
										<< "  --> pt = " << local_pT << std::endl
										<< "  --> pphi = " << local_pphi << std::endl
										<< "daughter_particle_id = " << daughter_particle_id << endl
										<< "parent_resonance_particle_id = " << parent_resonance_particle_id << endl
										<< "  --> Qfunc = " << Qfunc << endl
										<< "  --> n_body = " << n_body << endl
										<< "  --> gRES = " << gRES << endl
										<< "  --> Mres = " << Mres << endl
										<< "  --> mass = " << mass << endl
										<< "  --> Gamma = " << Gamma << endl
										<< "  --> br = " << br << endl
										<< "  --> m2 = " << m2 << endl
										<< "  --> m3 = " << m3 << endl << endl;
				exit(1);
			}
		}								// end of pT, pphi loops
	}										// end of nbody == 3



	//clean up
	Delete_resonance_running_sum_vectors();

	return;
}

inline void CorrelationFunction::set_to_zero(double * array, size_t arraylength)
{
	for (size_t arrayidx=0; arrayidx<arraylength; ++arrayidx) array[arrayidx] = 0.0;
}

//End of file
