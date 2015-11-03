#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>
#include<algorithm>

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

const double PTCHANGE = 1.0;

template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

CorrelationFunction::CorrelationFunction(particle_info* particle, particle_info* all_particles_in, int Nparticle_in,
					FO_surf* FOsurf_ptr, vector<int> chosen_resonances_in, int particle_idx, ofstream& myout)
{
	//set ofstream for output file
	global_out_stream_ptr = &myout;
	
	//particle information (both final-state particle used in HBT and all decay decay_channels)
	particle_name = particle->name;
	particle_mass = particle->mass;
	particle_sign = particle->sign;
	particle_gspin = particle->gspin;
	particle_id = particle_idx;
	target_particle_id = particle_id;
	particle_monval = particle->monval;
	S_prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;
	all_particles = all_particles_in;
	for (int icr = 0; icr < (int)chosen_resonances_in.size(); icr++)
		chosen_resonances.push_back(chosen_resonances_in[icr]);
	thermal_pions_only = false;
	Nparticle = Nparticle_in;
	read_in_all_dN_dypTdpTdphi = false;
	output_all_dN_dypTdpTdphi = true;
	currentfolderindex = -1;
	current_level_of_output = 0;

	Set_q_points();

	n_zeta_pts = 12;
	n_v_pts = 12;
	n_s_pts = 12;
	v_min = -1.;
	v_max = 1.;
	zeta_min = 0.;
	zeta_max = M_PI;
	
	//default: use delta_f in calculations
	use_delta_f = true;
	no_df_stem = "";
	current_FOsurf_ptr = FOsurf_ptr;

//****************************************************************************************************
//OLD CODE FOR READING IN SELECTED decay_channels...
	current_resonance_mass = 0.0;
	current_resonance_mu = 0.0;
	current_resonance_Gamma = 0.0;
	current_resonance_total_br = 0.0;
	current_resonance_decay_masses = new double [2];
	current_resonance_decay_masses[0] = 0.0;
	current_resonance_decay_masses[1] = 0.0;
	previous_resonance_particle_id = -1;
	previous_decay_channel_idx = -1;				//different for each decay channel
	previous_resonance_mass = 0.0;
	previous_resonance_Gamma = 0.0;
	previous_resonance_total_br = 0.0;
	if (chosen_resonances.size() == 0)
	{
		n_decay_channels = 1;
		n_resonance = 0;
		thermal_pions_only = true;
		if (VERBOSE > 0) *global_out_stream_ptr << "Thermal pion(+) only!" << endl;
		decay_channels = new decay_info [n_decay_channels];
	}
	else
	{
		//n_decay_channels is actually total number of decay channels which can generate pions
		//from chosen decay_channels
		n_decay_channels = get_number_of_decay_channels(chosen_resonances, all_particles);
		n_resonance = (int)chosen_resonances.size();
		if (VERBOSE > 0) *global_out_stream_ptr << "Computed n_decay_channels = " << n_decay_channels << endl
							<< "Computed n_resonance = " << n_resonance << endl;
		decay_channels = new decay_info [n_decay_channels];
		int temp_idx = 0;
		for (int icr = 0; icr < n_resonance; icr++)
		{
			particle_info particle_temp = all_particles[chosen_resonances[icr]];
			if (VERBOSE > 0) *global_out_stream_ptr << "Loading resonance: " << particle_temp.name
					<< ", chosen_resonances[" << icr << "] = " << chosen_resonances[icr] << endl;
			for (int idecay = 0; idecay < particle_temp.decays; idecay++)
			{
				if (VERBOSE > 0) *global_out_stream_ptr << "Current temp_idx = " << temp_idx << endl;
				if (temp_idx == n_decay_channels)	//i.e., all contributing decay channels have been loaded
					break;
				decay_channels[temp_idx].resonance_name = particle_temp.name;		// set name of resonance

				//check if effective branching is too small for inclusion in source variances
				bool effective_br_is_too_small = false;
				if (particle_temp.decays_effective_branchratio[idecay] <= 1.e-12)
					effective_br_is_too_small = true;

				decay_channels[temp_idx].resonance_particle_id = chosen_resonances[icr];	// set index of resonance in all_particles
				decay_channels[temp_idx].resonance_idx = icr;					// set index of resonance in chosen_resonances
				decay_channels[temp_idx].resonance_decay_masses = new double [Maxdecaypart];	// Maxdecaypart == 5
				decay_channels[temp_idx].resonance_decay_monvals = new double [Maxdecaypart];	// Maxdecaypart == 5
				decay_channels[temp_idx].resonance_decay_Gammas = new double [Maxdecaypart];	// Maxdecaypart == 5



				//*** SETTING RESONANCE DECAY MASSES DIFFERENTLY FOR NEW ANALYZE SF
				for (int ii = 0; ii < Maxdecaypart; ii++)
				{
					decay_channels[temp_idx].resonance_decay_monvals[ii] = particle_temp.decays_part[idecay][ii];
					if (particle_temp.decays_part[idecay][ii] == 0)
					{
						decay_channels[temp_idx].resonance_decay_masses[ii] = 0.0;
						decay_channels[temp_idx].resonance_decay_Gammas[ii] = 0.0;

					}
					else
					{
						int tempID = lookup_particle_id_from_monval(all_particles, Nparticle, particle_temp.decays_part[idecay][ii]);
						decay_channels[temp_idx].resonance_decay_masses[ii] = all_particles[tempID].mass;
						decay_channels[temp_idx].resonance_decay_Gammas[ii] = all_particles[tempID].width;

					}
				}
				decay_channels[temp_idx].resonance_mu = particle_temp.mu;
				decay_channels[temp_idx].resonance_gspin = particle_temp.gspin;
				decay_channels[temp_idx].resonance_sign = particle_temp.sign;
				decay_channels[temp_idx].resonance_mass = particle_temp.mass;
				decay_channels[temp_idx].nbody = abs(particle_temp.decays_Npart[idecay]);
				decay_channels[temp_idx].resonance_Gamma = particle_temp.width;
				decay_channels[temp_idx].resonance_total_br = particle_temp.decays_effective_branchratio[idecay];
				decay_channels[temp_idx].resonance_direct_br = particle_temp.decays_branchratio[idecay];

				
				//check if particle lifetime is too long for inclusion in source variances
				bool lifetime_is_too_long = false;
				//if (decay_channels[temp_idx].resonance_Gamma < hbarC / max_lifetime)
				//	lifetime_is_too_long = true;		//i.e., for lifetimes longer than 100 fm/c, skip decay channel

				if (VERBOSE > 0) *global_out_stream_ptr << "Resonance = " << decay_channels[temp_idx].resonance_name << ", decay channel " << idecay + 1
						<< ": mu=" << decay_channels[temp_idx].resonance_mu
						<< ", gs=" << decay_channels[temp_idx].resonance_gspin << ", sign=" << decay_channels[temp_idx].resonance_sign
						<< ", M=" << decay_channels[temp_idx].resonance_mass << ", nbody=" << decay_channels[temp_idx].nbody
						<< ", Gamma=" << decay_channels[temp_idx].resonance_Gamma << ", total br=" << decay_channels[temp_idx].resonance_total_br
						<< ", direct br=" << decay_channels[temp_idx].resonance_direct_br << endl;

				if (VERBOSE > 0) *global_out_stream_ptr << "Resonance = " << decay_channels[temp_idx].resonance_name << ": ";
				for (int decay_part_idx = 0; decay_part_idx < decay_channels[temp_idx].nbody; decay_part_idx++)
					if (VERBOSE > 0) *global_out_stream_ptr << "m[" << decay_part_idx << "] = "
						<< decay_channels[temp_idx].resonance_decay_masses[decay_part_idx] << "   "
						<< decay_channels[temp_idx].resonance_decay_monvals[decay_part_idx] << "   ";
				if (VERBOSE > 0) *global_out_stream_ptr << endl << endl;

				// if decay channel parent resonance is not too long-lived
				// and decay channel contains at least one target daughter particle,
				// include channel
				if (DO_ALL_DECAY_CHANNELS)
					decay_channels[temp_idx].include_channel = true;
				else
					decay_channels[temp_idx].include_channel = (!lifetime_is_too_long
																	&& !effective_br_is_too_small);

				temp_idx++;
			}
		}
	}
	
	//*****************************************************************
	// Only make dN_dypTdpTdphi_moments large enough to hold all necessary resonance, not decay channels
	//*****************************************************************
	dN_dypTdpTdphi_moments = new double ***** [Nparticle];
	ln_dN_dypTdpTdphi_moments = new double ***** [Nparticle];
	sign_of_dN_dypTdpTdphi_moments = new double ***** [Nparticle];
	for (int ir=0; ir<Nparticle; ir++)
	{
		dN_dypTdpTdphi_moments[ir] = new double **** [qnpts];
		ln_dN_dypTdpTdphi_moments[ir] = new double **** [qnpts];
		sign_of_dN_dypTdpTdphi_moments[ir] = new double **** [qnpts];
		for (int iq = 0; iq < qnpts; ++iq)
		{
			dN_dypTdpTdphi_moments[ir][iq] = new double *** [3];
			ln_dN_dypTdpTdphi_moments[ir][iq] = new double *** [3];
			sign_of_dN_dypTdpTdphi_moments[ir][iq] = new double *** [3];
			for (int iqax = 0; iqax < 3; ++iqax)
			{
				dN_dypTdpTdphi_moments[ir][iq][iqax] = new double ** [2];
				ln_dN_dypTdpTdphi_moments[ir][iq][iqax] = new double ** [2];
		        sign_of_dN_dypTdpTdphi_moments[ir][iq][iqax] = new double ** [2];
				for (int itrig = 0; itrig < 2; ++itrig)
				{
					dN_dypTdpTdphi_moments[ir][iq][iqax][itrig] = new double * [n_interp_pT_pts];
					ln_dN_dypTdpTdphi_moments[ir][iq][iqax][itrig] = new double * [n_interp_pT_pts];
					sign_of_dN_dypTdpTdphi_moments[ir][iq][iqax][itrig] = new double * [n_interp_pT_pts];
					for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
					{
						dN_dypTdpTdphi_moments[ir][iq][iqax][itrig][ipT] = new double [n_interp_pphi_pts];
                        ln_dN_dypTdpTdphi_moments[ir][iq][iqax][itrig][ipT] = new double [n_interp_pphi_pts];
		                sign_of_dN_dypTdpTdphi_moments[ir][iq][iqax][itrig][ipT] = new double [n_interp_pphi_pts];
                        for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
						{
							dN_dypTdpTdphi_moments[ir][iq][iqax][itrig][ipT][ipphi] = 0.0;
                            ln_dN_dypTdpTdphi_moments[ir][iq][iqax][itrig][ipT][ipphi] = 0.0;
		                    sign_of_dN_dypTdpTdphi_moments[ir][iq][iqax][itrig][ipT][ipphi] = 0.0;
                        }
					}
				}
			}
		}
	}

	// initialize spectra and correlation function arrays
	spectra = new double ** [Nparticle];
	for (int ir=0; ir<Nparticle; ir++)
	{
		spectra[ir] = new double * [n_interp_pT_pts];
		for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
		{
			spectra[ir][ipT] = new double [n_interp_pphi_pts];
			for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
				spectra[ir][ipT][ipphi] = 0.0;
		}
	}

	CFvals = new double *** [n_interp_pT_pts];
	for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
	{
		CFvals[ipT] = new double ** [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			CFvals[ipT][ipphi] = new double * [qnpts];
			for (int iq = 0; iq < qnpts; ++iq)
			{
				CFvals[ipT][ipphi][iq] = new double [3];
				for (int iqax = 0; iqax < 3; ++iqax)
					CFvals[ipT][ipphi][iq][iqax] = 0.0;
			}
		}
	}

	// set-up integration points for resonance integrals
	s_pts = new double * [n_decay_channels];
	s_wts = new double * [n_decay_channels];
	v_pts = new double [n_v_pts];
	v_wts = new double [n_v_pts];
	zeta_pts = new double [n_zeta_pts];
	zeta_wts = new double [n_zeta_pts];
	for (int idc=0; idc<n_decay_channels; idc++)
	{
		s_pts[idc] = new double [n_s_pts];
		s_wts[idc] = new double [n_s_pts];
	}

	//initialize all gaussian points for resonance integrals
	//syntax: int gauss_quadrature(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[])
	gauss_quadrature(n_zeta_pts, 1, 0.0, 0.0, zeta_min, zeta_max, zeta_pts, zeta_wts);
	gauss_quadrature(n_v_pts, 1, 0.0, 0.0, v_min, v_max, v_pts, v_wts);
	for (int idc = 0; idc < n_decay_channels; idc++)
	{
		//cerr << "working on resonance #" << idc << "..." << endl;
		double Gamma_temp = decay_channels[idc].resonance_Gamma;
		double m2_temp = decay_channels[idc].resonance_decay_masses[0];
		double m3_temp = decay_channels[idc].resonance_decay_masses[1];
		double M_temp = decay_channels[idc].resonance_mass;
		double s_min_temp = (m2_temp + m3_temp)*(m2_temp + m3_temp);
		double s_max_temp = (M_temp - particle_mass)*(M_temp - particle_mass);
		// N.B. - this is only really necessary for 3-body decays,
		//			but doesn't cause any problems for 2-body and is easier/simpler to code...
		gauss_quadrature(n_s_pts, 1, 0.0, 0.0, s_min_temp, s_max_temp, s_pts[idc], s_wts[idc]);
	}	

   //single particle spectra for plane angle determination
   SP_pT = new double [n_SP_pT];
   SP_pT_weight = new double [n_SP_pT];
   gauss_quadrature(n_SP_pT, 1, 0.0, 0.0, SP_pT_min, SP_pT_max, SP_pT, SP_pT_weight);
   SP_pphi = new double [n_SP_pphi];
   SP_pphi_weight = new double [n_SP_pphi];
   gauss_quadrature(n_SP_pphi, 1, 0.0, 0.0, 0.0, 2.*M_PI, SP_pphi, SP_pphi_weight);
   SP_p_y = 0.0e0;

//initialize and set evenly spaced grid of px-py points in transverse plane,
//and corresponding p0 and pz points
	SPinterp_pT = new double [n_interp_pT_pts];
	SPinterp_pphi = new double [n_interp_pphi_pts];
	sin_SPinterp_pphi = new double [n_interp_pphi_pts];
	cos_SPinterp_pphi = new double [n_interp_pphi_pts];
	SPinterp_p0 = new double* [n_interp_pT_pts];
	SPinterp_pz = new double* [n_interp_pT_pts];
	double * dummywts3 = new double [n_interp_pT_pts];
	double * dummywts4 = new double [n_interp_pphi_pts];
	for(int ipt=0; ipt<n_interp_pT_pts; ipt++)
	{
		SPinterp_p0[ipt] = new double [eta_s_npts];
		SPinterp_pz[ipt] = new double [eta_s_npts];
	}
	if (UNIFORM_SPACING)
	{
		//use uniformly spaced points in transverse momentum to make
		//interpolation simpler/faster
		for(int ipt=0; ipt<n_interp_pT_pts; ipt++)
			SPinterp_pT[ipt] = interp_pT_min + (double)ipt*Del2_pT;
		for(int ipphi=0; ipphi<n_interp_pphi_pts; ipphi++)
		{
			SPinterp_pphi[ipphi] = interp_pphi_min + (double)ipphi*Del2_pphi;
			sin_SPinterp_pphi[ipphi] = sin(SPinterp_pphi[ipphi]);
			cos_SPinterp_pphi[ipphi] = cos(SPinterp_pphi[ipphi]);
		}
	}
	else
	{
		//try just gaussian points...
		//syntax:
		//int gauss_quadrature(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[]) 
		gauss_quadrature(n_interp_pT_pts, 5, 0.0, 0.0, 0.0, 13.0, SPinterp_pT, dummywts3);	//use this one to agree with iS.e
		gauss_quadrature(n_interp_pphi_pts, 1, 0.0, 0.0, interp_pphi_min, interp_pphi_max, SPinterp_pphi, dummywts4);
		for(int ipphi=0; ipphi<n_interp_pphi_pts; ipphi++)
		{
			//SPinterp_pphi[ipphi] = interp_pphi_min + (double)ipphi*Del2_pphi;
			sin_SPinterp_pphi[ipphi] = sin(SPinterp_pphi[ipphi]);
			cos_SPinterp_pphi[ipphi] = cos(SPinterp_pphi[ipphi]);
		}
	}

	dN_dypTdpTdphi = new double* [n_SP_pT];
	cosine_iorder = new double* [n_SP_pT];
	sine_iorder = new double* [n_SP_pT];
	for(int i=0; i<n_SP_pT; i++)
	{
		dN_dypTdpTdphi[i] = new double [n_SP_pphi];
		cosine_iorder[i] = new double [n_order];
		sine_iorder[i] = new double [n_order];
	}
	dN_dydphi = new double [n_SP_pphi];
	dN_dypTdpT = new double [n_SP_pT];
	pTdN_dydphi = new double [n_SP_pphi];
	for(int i=0; i<n_SP_pphi; i++)
	{
		dN_dydphi[i] = 0.0e0;
		pTdN_dydphi[i] = 0.0e0;
		for(int j=0; j<n_SP_pT; j++)
			dN_dypTdpTdphi[j][i] = 0.0e0;
	}
	for(int i=0; i<n_SP_pT; i++)
	{
		dN_dypTdpT[i] = 0.0e0;
		for(int j=0; j<n_order; j++)
		{
			cosine_iorder[i][j] = 0.0e0;
			sine_iorder[i][j] = 0.0e0;
		}
	}
	plane_angle = new double [n_order];

	//pair momentum
	K_T = new double [n_localp_T];
	double dK_T = (localp_T_max - localp_T_min)/(n_localp_T - 1 + 1e-100);
	for(int i=0; i<n_localp_T; i++) K_T[i] = localp_T_min + i*dK_T;
	//K_y = p_y;
	K_y = 0.;
	ch_K_y = cosh(K_y);
	sh_K_y = sinh(K_y);
	beta_l = sh_K_y/ch_K_y;
	K_phi = new double [n_localp_phi];
	K_phi_weight = new double [n_localp_phi];
	gauss_quadrature(n_localp_phi, 1, 0.0, 0.0, localp_phi_min, localp_phi_max, K_phi, K_phi_weight);

	//spatial rapidity grid
	eta_s = new double [eta_s_npts];
	eta_s_weight = new double [eta_s_npts];
	gauss_quadrature(eta_s_npts, 1, 0.0, 0.0, eta_s_i, eta_s_f, eta_s, eta_s_weight);
	ch_eta_s = new double [eta_s_npts];
	sh_eta_s = new double [eta_s_npts];
	for (int ieta = 0; ieta < eta_s_npts; ieta++)
	{
		ch_eta_s[ieta] = cosh(eta_s[ieta]);
		sh_eta_s[ieta] = sinh(eta_s[ieta]);
	}

	R2_side = new double* [n_localp_T];
	R2_side_C = new double* [n_localp_T];
	R2_side_S = new double* [n_localp_T];
	R2_out = new double* [n_localp_T];
	R2_out_C = new double* [n_localp_T];
	R2_out_S = new double* [n_localp_T];
	R2_long = new double* [n_localp_T];
	R2_long_C = new double* [n_localp_T];
	R2_long_S = new double* [n_localp_T];
	R2_outside = new double* [n_localp_T];
	R2_outside_C = new double* [n_localp_T];
	R2_outside_S = new double* [n_localp_T];
	R2_sidelong = new double* [n_localp_T];
	R2_sidelong_C = new double* [n_localp_T];
	R2_sidelong_S = new double* [n_localp_T];
	R2_outlong = new double* [n_localp_T];
	R2_outlong_C = new double* [n_localp_T];
	R2_outlong_S = new double* [n_localp_T];

	R2_side_err = new double* [n_localp_T];
	R2_out_err = new double* [n_localp_T];
	R2_long_err = new double* [n_localp_T];
	R2_outside_err = new double* [n_localp_T];
	R2_sidelong_err = new double* [n_localp_T];
	R2_outlong_err = new double* [n_localp_T];

	lambda_Correl = new double * [n_localp_T];
	lambda_Correl_err = new double * [n_localp_T];

	for(int i=0; i<n_localp_T; i++)
	{
		R2_side[i] = new double [n_localp_phi];
		R2_side_C[i] = new double [n_order];
		R2_side_S[i] = new double [n_order];
		R2_out[i] = new double [n_localp_phi];
		R2_out_C[i] = new double [n_order];
		R2_out_S[i] = new double [n_order];
		R2_outside[i] = new double [n_localp_phi];
		R2_outside_C[i] = new double [n_order];
		R2_outside_S[i] = new double [n_order];
		R2_long[i] = new double [n_localp_phi];
		R2_long_C[i] = new double [n_order];
		R2_long_S[i] = new double [n_order];
		R2_sidelong[i] = new double [n_localp_phi];
		R2_sidelong_C[i] = new double [n_order];
		R2_sidelong_S[i] = new double [n_order];
		R2_outlong[i] = new double [n_localp_phi];
		R2_outlong_C[i] = new double [n_order];
		R2_outlong_S[i] = new double [n_order];

		R2_side_err[i] = new double [n_localp_phi];
		R2_out_err[i] = new double [n_localp_phi];
		R2_long_err[i] = new double [n_localp_phi];
		R2_outside_err[i] = new double [n_localp_phi];
		R2_sidelong_err[i] = new double [n_localp_phi];
		R2_outlong_err[i] = new double [n_localp_phi];

		lambda_Correl[i] = new double [n_localp_phi];
		lambda_Correl_err[i] = new double [n_localp_phi];
	}

	//initialize all source variances and HBT radii/coeffs
	for(int i=0; i<n_localp_T; i++)
	{
		for(int j=0; j<n_localp_phi; j++)
		{
			R2_side[i][j] = 0.;
			R2_out[i][j] = 0.;
			R2_long[i][j] = 0.;
			R2_outside[i][j] = 0.;
			R2_sidelong[i][j] = 0.;
			R2_outlong[i][j] = 0.;

			R2_side_err[i][j] = 0.;
			R2_out_err[i][j] = 0.;
			R2_long_err[i][j] = 0.;
			R2_outside_err[i][j] = 0.;
			R2_sidelong_err[i][j] = 0.;
			R2_outlong_err[i][j] = 0.;

			lambda_Correl[i][j] = 0.0;
			lambda_Correl_err[i][j] = 0.0;
		}
		for(int j=0; j<n_order; j++)
		{
			R2_side_C[i][j] = 0.;
			R2_side_S[i][j] = 0.;
			R2_out_C[i][j] = 0.;
			R2_out_S[i][j] = 0.;
			R2_outside_C[i][j] = 0.;
			R2_outside_S[i][j] = 0.;
			R2_long_C[i][j] = 0.;
			R2_long_S[i][j] = 0.;
			R2_sidelong_C[i][j] = 0.;
			R2_sidelong_S[i][j] = 0.;
			R2_outlong_C[i][j] = 0.;
			R2_outlong_S[i][j] = 0.;
		}
	}

   return;
}

void CorrelationFunction::Update_sourcefunction(particle_info* particle, int FOarray_length, int particle_idx)
{
   //particle information
   particle_name = particle->name;
   particle_mass = particle->mass;
   particle_sign = particle->sign;
   particle_gspin = particle->gspin;
   particle_id = particle_idx;

	*global_out_stream_ptr << "Inside Update_sourcefunction(...): using fraction_of_resonances = " << fraction_of_resonances << endl;

   //erase contents of single - and two-particle spectra
   for(int i=0; i<n_SP_pphi; i++)
   {
      dN_dydphi[i] = 0.0e0;
      pTdN_dydphi[i] = 0.0e0;
      for(int j=0; j<n_SP_pT; j++) dN_dypTdpTdphi[j][i] = 0.0e0;
   }
   //erase anisotropic flows
   for(int i=0; i<n_SP_pT; i++)
   for(int j=0; j<n_order; j++)
   {
      cosine_iorder[i][j] = 0.0e0;
      sine_iorder[i][j] = 0.0e0;
   }

   FO_length = FOarray_length;

//reset only EBE source variances and EBE HBT radii/coeffs
for(int i=0; i<n_localp_T; i++)
{
	for(int j=0; j<n_order; j++)
	{
		R2_side_C[i][j] = 0.;
		R2_side_S[i][j] = 0.;
		R2_out_C[i][j] = 0.;
		R2_out_S[i][j] = 0.;
		R2_outside_C[i][j] = 0.;
		R2_outside_S[i][j] = 0.;
		R2_long_C[i][j] = 0.;
		R2_long_S[i][j] = 0.;
		R2_sidelong_C[i][j] = 0.;
		R2_sidelong_S[i][j] = 0.;
		R2_outlong_C[i][j] = 0.;
		R2_outlong_S[i][j] = 0.;
	}
}

   return;
}

CorrelationFunction::~CorrelationFunction()
{
   delete[] SP_pT;
   delete[] SP_pT_weight;
   delete[] SP_pphi;
   delete[] SP_pphi_weight;
   delete[] dN_dydphi;
   delete[] dN_dypTdpT;
   delete[] pTdN_dydphi;
   for(int i=0; i<n_SP_pT; i++)
   {
      delete[] dN_dypTdpTdphi[i];
      delete[] cosine_iorder[i];
      delete[] sine_iorder[i];
   }
   delete[] dN_dypTdpTdphi;
   delete[] cosine_iorder;
   delete[] sine_iorder;
   delete[] plane_angle;

   delete[] K_T;
   delete[] K_phi;
   delete[] K_phi_weight;
   delete[] eta_s;
   delete[] eta_s_weight;

   for(int i=0; i<n_localp_T; i++)
   {
      delete[] R2_side[i];
   }

   delete[] R2_side;

   return;
}

void CorrelationFunction::Allocate_decay_channel_info()
{
	if (VERBOSE > 2) *global_out_stream_ptr << "Reallocating memory for decay channel information..." << endl;
	VEC_n2_v_factor = new double [n_v_pts];
	VEC_n2_zeta_factor = new double * [n_v_pts];
	VEC_n2_P_Y = new double [n_v_pts];
	VEC_n2_MTbar = new double [n_v_pts];
	VEC_n2_DeltaMT = new double [n_v_pts];
	VEC_n2_MTp = new double [n_v_pts];
	VEC_n2_MTm = new double [n_v_pts];
	VEC_n2_MT = new double * [n_v_pts];
	VEC_n2_PPhi_tilde = new double * [n_v_pts];
	VEC_n2_PPhi_tildeFLIP = new double * [n_v_pts];
	VEC_n2_PT = new double * [n_v_pts];
	for(int iv = 0; iv < n_v_pts; ++iv)
	{
		VEC_n2_MT[iv] = new double [n_zeta_pts];
		VEC_n2_PPhi_tilde[iv] = new double [n_zeta_pts];
		VEC_n2_PPhi_tildeFLIP[iv] = new double [n_zeta_pts];
		VEC_n2_PT[iv] = new double [n_zeta_pts];
		VEC_n2_zeta_factor[iv] = new double [n_zeta_pts];
	}
	NEW_s_pts = new double [n_s_pts];
	NEW_s_wts = new double [n_s_pts];
	VEC_pstar = new double [n_s_pts];
	VEC_Estar = new double [n_s_pts];
	VEC_DeltaY = new double [n_s_pts];
	VEC_g_s = new double [n_s_pts];
	VEC_s_factor = new double [n_s_pts];
	VEC_v_factor = new double * [n_s_pts];
	VEC_zeta_factor = new double ** [n_s_pts];
	VEC_Yp = new double [n_s_pts];
	VEC_Ym = new double [n_s_pts];
	VEC_P_Y = new double * [n_s_pts];
	VEC_MTbar = new double * [n_s_pts];
	VEC_DeltaMT = new double * [n_s_pts];
	VEC_MTp = new double * [n_s_pts];
	VEC_MTm = new double * [n_s_pts];
	VEC_MT = new double ** [n_s_pts];
	VEC_PPhi_tilde = new double ** [n_s_pts];
	VEC_PPhi_tildeFLIP = new double ** [n_s_pts];
	VEC_PT = new double ** [n_s_pts];
	for(int is = 0; is < n_s_pts; ++is)
	{
		VEC_v_factor[is] = new double [n_v_pts];
		VEC_zeta_factor[is] = new double * [n_v_pts];
		VEC_P_Y[is] = new double [n_v_pts];
		VEC_MTbar[is] = new double [n_v_pts];
		VEC_DeltaMT[is] = new double [n_v_pts];
		VEC_MTp[is] = new double [n_v_pts];
		VEC_MTm[is] = new double [n_v_pts];
		VEC_MT[is] = new double * [n_v_pts];
		VEC_PPhi_tilde[is] = new double * [n_v_pts];
		VEC_PPhi_tildeFLIP[is] = new double * [n_v_pts];
		VEC_PT[is] = new double * [n_v_pts];
		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			VEC_MT[is][iv] = new double [n_zeta_pts];
			VEC_PPhi_tilde[is][iv] = new double [n_zeta_pts];
			VEC_PPhi_tildeFLIP[is][iv] = new double [n_zeta_pts];
			VEC_PT[is][iv] = new double [n_zeta_pts];
			VEC_zeta_factor[is][iv] = new double [n_zeta_pts];
		}
	}
	if (VERBOSE > 2) *global_out_stream_ptr << "Reallocated memory for decay channel information." << endl;

	return;
}

void CorrelationFunction::Delete_decay_channel_info()
{
	if (VERBOSE > 2) *global_out_stream_ptr << "Deleting memory for decay channel information..." << endl;
	for(int iv = 0; iv < n_v_pts; ++iv)
	{
		delete [] VEC_n2_MT[iv];
		delete [] VEC_n2_PPhi_tilde[iv];
		delete [] VEC_n2_PPhi_tildeFLIP[iv];
		delete [] VEC_n2_PT[iv];
		delete [] VEC_n2_zeta_factor[iv];
	}
	delete [] VEC_n2_v_factor;
	delete [] VEC_n2_zeta_factor;
	delete [] VEC_n2_P_Y;
	delete [] VEC_n2_MTbar ;
	delete [] VEC_n2_DeltaMT;
	delete [] VEC_n2_MTp;
	delete [] VEC_n2_MTm;
	delete [] VEC_n2_MT;
	delete [] VEC_n2_PPhi_tilde;
	delete [] VEC_n2_PPhi_tildeFLIP;
	delete [] VEC_n2_PT;

	for(int is = 0; is < n_s_pts; ++is)
	{
		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			delete [] VEC_MT[is][iv];
			delete [] VEC_PPhi_tilde[is][iv];
			delete [] VEC_PPhi_tildeFLIP[is][iv];
			delete [] VEC_PT[is][iv];
			delete [] VEC_zeta_factor[is][iv];
		}
		delete [] VEC_v_factor[is];
		delete [] VEC_zeta_factor[is];
		delete [] VEC_P_Y[is];
		delete [] VEC_MTbar[is];
		delete [] VEC_DeltaMT[is];
		delete [] VEC_MTp[is];
		delete [] VEC_MTm[is];
		delete [] VEC_MT[is];
		delete [] VEC_PPhi_tilde[is];
		delete [] VEC_PPhi_tildeFLIP[is];
		delete [] VEC_PT[is];
	}
	delete [] NEW_s_pts;
	delete [] NEW_s_wts;
	delete [] VEC_pstar;
	delete [] VEC_Estar;
	delete [] VEC_DeltaY;
	delete [] VEC_g_s;
	delete [] VEC_s_factor;
	delete [] VEC_v_factor;
	delete [] VEC_zeta_factor;
	delete [] VEC_Yp;
	delete [] VEC_Ym;
	delete [] VEC_P_Y;
	delete [] VEC_MTbar;
	delete [] VEC_DeltaMT;
	delete [] VEC_MTp;
	delete [] VEC_MTm;
	delete [] VEC_MT;
	delete [] VEC_PPhi_tilde;
	delete [] VEC_PPhi_tildeFLIP;
	delete [] VEC_PT;
	if (VERBOSE > 2) *global_out_stream_ptr << "Deleted memory for decay channel information." << endl;

	return;
}

void CorrelationFunction::Set_q_points()
{
	q_pts = new double [qnpts];
	for (int iq = 0; iq < qnpts; ++iq)
		q_pts[iq] = init_q + (double)iq * delta_q;
	q_axes = new double [3];
}

bool CorrelationFunction::fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

//print output to output filestream, one line at a time
void CorrelationFunction::Set_ofstream(ofstream& myout)
{
	global_out_stream_ptr = &myout;

	return;
}

//print output to output filestream, one line at a time
void CorrelationFunction::Set_path(string localpath)
{
	global_path = localpath;

	return;
}

void CorrelationFunction::Set_runfolder(string localrunfolder)
{
	global_runfolder = localrunfolder;

	return;
}

void CorrelationFunction::Set_resultsfolder_stem(string usrdef_resultsfolder_stem)
{
	global_resultsfolder_stem = usrdef_resultsfolder_stem;

	return;
}

void CorrelationFunction::Set_use_delta_f(bool usrdef_usedeltaf)
{
	use_delta_f = usrdef_usedeltaf;
	if (!use_delta_f)
		no_df_stem = "_no_df";
	return;
}

void CorrelationFunction::Set_particle_mass(double usrdef_particle_mass)
{
	particle_mass = usrdef_particle_mass;
	return;
}

void CorrelationFunction::Set_current_FOsurf_ptr(FO_surf* FOsurf_ptr)
{
	current_FOsurf_ptr = FOsurf_ptr;
	return;
}

void CorrelationFunction::Get_current_decay_string(int dc_idx, string * decay_string)
{
	// N.B. - dc_idx == 0 is thermal pion(+)s in calling loop, dc_idx > 0 gives resonance decays
	//      ==> need to use dc_idx - 1 here
	*decay_string = decay_channels[dc_idx - 1].resonance_name + " --->> ";
	int temp_monval, tempID;
	for (int decay_part_idx = 0; decay_part_idx < decay_channels[dc_idx - 1].nbody; decay_part_idx++)
	{
		temp_monval = decay_channels[dc_idx - 1].resonance_decay_monvals[decay_part_idx];
		//if (VERBOSE > 0) *global_out_stream_ptr << "Get_current_decay_string(): temp_monval = " << temp_monval << endl;
		if (temp_monval == 0)
			continue;
		else
		{
			tempID = lookup_particle_id_from_monval(all_particles, Nparticle, temp_monval);
			*decay_string += all_particles[tempID].name;
			if (decay_part_idx < decay_channels[dc_idx - 1].nbody - 1) *decay_string += " + ";
		}
	}
	return;
}

int CorrelationFunction::lookup_resonance_idx_from_particle_id(int pid)
{
	// pid - particle index in all_particles array
	// looks up location in chosen_resonances of given value particle_id
	int result = -1;

	for (int ii = 0; ii < (int)chosen_resonances.size(); ii++)
	{
		if (chosen_resonances[ii] == pid)
		{
			result = ii;
			break;
		}
	}

	// if pid is not one of the chosen_resonances, is not the target daughter (pion(+)), is not stable and has a non-zero effective branching ratio
	if (result < 0 && pid != particle_id && all_particles[pid].stable == 0 && all_particles[pid].effective_branchratio >= 1.e-12)
	{
		*global_out_stream_ptr << " *** lookup_resonance_idx_from_particle_id(): Particle_id = " << pid
					<< " (" << all_particles[pid].name <<") not found in chosen_resonances!" << endl
					<< " *** br = " << all_particles[pid].effective_branchratio << endl;
	}
	return (result);
}

void CorrelationFunction::Allocate_resonance_running_sum_vectors()
{
    ssum_vec = new double ** [qnpts];
    vsum_vec = new double ** [qnpts];
    zetasum_vec = new double ** [qnpts];
    Csum_vec = new double ** [qnpts];
	for (int iq = 0; iq < qnpts; ++iq)
	{
        ssum_vec[iq] = new double * [3];
        vsum_vec[iq] = new double * [3];
        zetasum_vec[iq] = new double * [3];
        Csum_vec[iq] = new double * [3];
	    for (int iqax = 0; iqax < 3; ++iqax)
	    {
            ssum_vec[iq][iqax] = new double [2];
            vsum_vec[iq][iqax] = new double [2];
            zetasum_vec[iq][iqax] = new double [2];
            Csum_vec[iq][iqax] = new double [2];
	        for (int itrig = 0; itrig < 2; ++itrig)
	        {
		        ssum_vec[iq][iqax][itrig] = 0.0;
		        vsum_vec[iq][iqax][itrig] = 0.0;
		        zetasum_vec[iq][iqax][itrig] = 0.0;
		        Csum_vec[iq][iqax][itrig] = 0.0;
            }
        }
    }
}

void CorrelationFunction::Delete_resonance_running_sum_vectors()
{
	for (int iq = 0; iq < qnpts; ++iq)
	{
	    for (int iqax = 0; iqax < 3; ++iqax)
	    {
            delete [] ssum_vec[iq][iqax];
            delete [] vsum_vec[iq][iqax];
            delete [] zetasum_vec[iq][iqax];
            delete [] Csum_vec[iq][iqax];
        }
        delete [] ssum_vec[iq];
        delete [] vsum_vec[iq];
        delete [] zetasum_vec[iq];
        delete [] Csum_vec[iq];
    }
    delete [] ssum_vec;
    delete [] vsum_vec;
    delete [] zetasum_vec;
    delete [] Csum_vec;
}

void CorrelationFunction::Zero_resonance_running_sum_vector(double *** vec)
{
	for (int iq = 0; iq < qnpts; ++iq)
	for (int iqax = 0; iqax < 3; ++iqax)
	for (int itrig = 0; itrig < 2; ++itrig)
		vec[iq][iqax][itrig] = 0.0;
}

inline double CorrelationFunction::lin_int(double x_m_x1, double one_by_x2_m_x1, double f1, double f2)
{
	return ( f1 + (f2 - f1) * x_m_x1 * one_by_x2_m_x1 );
}

void CorrelationFunction::Edndp3(double ptr, double phir, double *** results)
{
	double phi0, phi1;
	double f1, f2;

	int npphi_max = n_interp_pphi_pts - 1;
	int npT_max = n_interp_pT_pts - 1;

	// locate pT interval
	int npt = 1;
	while ((ptr > SPinterp_pT[npt]) &&
			(npt < npT_max)) ++npt;
	double pT0 = SPinterp_pT[npt-1];
	double pT1 = SPinterp_pT[npt];

	// locate pphi interval
	int nphi = 1, nphim1 = 0;
	if(phir < SPinterp_pphi[0])			//if angle is less than minimum angle grid point
	{
		phi0 = SPinterp_pphi[npphi_max] - 2. * M_PI;
		phi1 = SPinterp_pphi[0];
		nphi = 0;
		nphim1 = npphi_max;
	}
	else if(phir > SPinterp_pphi[npphi_max])	//if angle is greater than maximum angle grid point
	{
		phi0 = SPinterp_pphi[npphi_max];
		phi1 = SPinterp_pphi[0] + 2. * M_PI;
		nphi = 0;
		nphim1 = npphi_max;
	}
	else						//if angle is within grid range
	{
		while ((phir > SPinterp_pphi[nphi]) &&
				(nphi < npphi_max)) ++nphi;
		nphim1 = nphi - 1;
		phi0 = SPinterp_pphi[nphim1];
		phi1 = SPinterp_pphi[nphi];
	}

	if (pT0==pT1 || phi0==phi1)
	{
		cerr << "ERROR in Edndp3(): pT and/or pphi values equal!" << endl;
		exit(1);
	}

	double one_by_pTdiff = 1./(pT1 - pT0), one_by_pphidiff = 1./(phi1 - phi0);

	for (int iq = 0; iq < qnpts; ++iq)
	{
		double **** temp_res_sign_info = res_sign_info[iq];
		double **** temp_res_log_info = res_log_info[iq];
		double **** temp_res_moments_info = res_moments_info[iq];
		for (int iqax = 0; iqax < 3; ++iqax)
		{
			double *** temp2_res_sign_info = temp_res_sign_info[iqax];
			double *** temp2_res_log_info = temp_res_log_info[iqax];
			double *** temp2_res_moments_info = temp_res_moments_info[iqax];
			for (int itrig = 0; itrig < 2; ++itrig)
			{
				double ** temp3_res_sign_info = temp2_res_sign_info[itrig];
				double ** temp3_res_log_info = temp2_res_log_info[itrig];
				double ** temp3_res_moments_info = temp2_res_moments_info[itrig];
				// interpolate over pT values first
				if(ptr > PTCHANGE)				// if pT interpolation point is larger than PTCHANGE (currently 1.0 GeV)
				{
					double sign_of_f11 = temp3_res_sign_info[npt-1][nphim1];
					double sign_of_f12 = temp3_res_sign_info[npt-1][nphi];
					double sign_of_f21 = temp3_res_sign_info[npt][nphim1];
					double sign_of_f22 = temp3_res_sign_info[npt][nphi];
			
					//*******************************************************************************************************************
					// set f1 first
					//*******************************************************************************************************************
					// if using extrapolation and spectra at pT1 has larger magnitude than at pT0, just return zero
					if (ptr > pT1 && temp3_res_log_info[npt][nphim1] > temp3_res_log_info[npt-1][nphim1])
						f1 = 0.0;
					else if (sign_of_f11 * sign_of_f21 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
						f1 = sign_of_f11 * exp( lin_int(ptr-pT0, one_by_pTdiff, temp3_res_log_info[npt-1][nphim1], temp3_res_log_info[npt][nphim1]) );
					else					// otherwise, just interpolate original vals
						f1 = lin_int(ptr-pT0, one_by_pTdiff, temp3_res_moments_info[npt-1][nphim1], temp3_res_moments_info[npt][nphim1]);
		
					//*******************************************************************************************************************
					// set f2 next
					//*******************************************************************************************************************
					if (ptr > pT1 && temp3_res_log_info[npt][nphi] > temp3_res_log_info[npt-1][nphi])
						f2 = 0.0;
					else if (sign_of_f12 * sign_of_f22 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
						f2 = sign_of_f12 * exp( lin_int(ptr-pT0, one_by_pTdiff, temp3_res_log_info[npt-1][nphi], temp3_res_log_info[npt][nphi]) );
					else					// otherwise, just interpolate original vals
						f2 = lin_int(ptr-pT0, one_by_pTdiff, temp3_res_moments_info[npt-1][nphi], temp3_res_moments_info[npt][nphi]);
					//*******************************************************************************************************************
				}
				else						// if pT is smaller than PTCHANGE, just use linear interpolation, no matter what
				{
					f1 = lin_int(ptr-pT0, one_by_pTdiff, temp3_res_moments_info[npt-1][nphim1], temp3_res_moments_info[npt][nphim1]);
					f2 = lin_int(ptr-pT0, one_by_pTdiff, temp3_res_moments_info[npt-1][nphi], temp3_res_moments_info[npt][nphi]);
				}
			
				// now, interpolate f1 and f2 over the pphi direction
				results[iq][iqax][itrig] += lin_int(phir-phi0, one_by_pphidiff, f1, f2);
			
				if ( isnan( results[iq][iqax][itrig] ) )
				{
					*global_out_stream_ptr << "ERROR in Edndp3(double, double, double*): problems encountered!" << endl
						<< "results[" << iq << "][" << iqax << "][" << itrig << "] = " << setw(8) << setprecision(15) << results[iq][iqax][itrig] << endl
						<< "  --> ptr = " << ptr << endl
						<< "  --> pt0 = " << pT0 << endl
						<< "  --> pt1 = " << pT1 << endl
						<< "  --> phir = " << phir << endl
						<< "  --> phi0 = " << phi0 << endl
						<< "  --> phi1 = " << phi1 << endl
						<< "  --> f11 = " << temp3_res_moments_info[npt-1][nphim1] << endl
						<< "  --> f12 = " << temp3_res_moments_info[npt-1][nphi] << endl
						<< "  --> f21 = " << temp3_res_moments_info[npt][nphim1] << endl
						<< "  --> f22 = " << temp3_res_moments_info[npt][nphi] << endl
						<< "  --> f1 = " << f1 << endl
						<< "  --> f2 = " << f2 << endl;
					exit(1);
				}
			}
		}
	}

	return;
}

//End of file
