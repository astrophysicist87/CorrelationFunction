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

void replace_parentheses(std::string & tempstring)
{
	int len = tempstring.length();
	for(unsigned int i = 0; i < len; i++)
	{
		char c = tempstring[i];
		if (c == '(' || c == ')')
			tempstring[i] = '_';
	}
	
	if (tempstring[len - 1] == '_')
		tempstring.erase( len - 1 );
	
	return;
}

void CorrelationFunction::Output_results(int folderindex)
{
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputHBT;
	outputHBT.open(filename_stream_HBT.str().c_str());
	ostringstream filename_stream_HBTcfs;
	filename_stream_HBTcfs << global_path << "/HBTradii_cfs_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputHBTcoeffs(filename_stream_HBTcfs.str().c_str());

	for(int iKT = 0; iKT < n_localp_T; iKT++)
	{
		for(int Morder=0; Morder<n_order; Morder++)
		{
			outputHBTcoeffs << folderindex << "  " << K_T[iKT] << "  " << Morder
				<< "  " << R2_side_C[iKT][Morder] << "   " << R2_side_S[iKT][Morder] << "  " << R2_out_C[iKT][Morder] << "  " << R2_out_S[iKT][Morder]
				<< "  " << R2_outside_C[iKT][Morder] << "   " << R2_outside_S[iKT][Morder] << "  " << R2_long_C[iKT][Morder] << "  " << R2_long_S[iKT][Morder]
				<< "  " << R2_sidelong_C[iKT][Morder] << "   " << R2_sidelong_S[iKT][Morder] << "  " << R2_outlong_C[iKT][Morder] << "  " << R2_outlong_S[iKT][Morder] << endl;
		}
		for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
		{
			outputHBT << folderindex << "  " << K_T[iKT] << "  " << K_phi[iKphi]
				<< "  " << R2_side[iKT][iKphi] << "  " << R2_out[iKT][iKphi]
				<< "  " << R2_outside[iKT][iKphi] << "  " << R2_long[iKT][iKphi]
				<< "  " << R2_sidelong[iKT][iKphi] << "  " << R2_outlong[iKT][iKphi] << endl;
		}
	}

	outputHBT.close();

	return;
}

void CorrelationFunction::Output_Correlationfunction_1D(int folderindex)
{
	ostringstream oCorrFunc_1D_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);
	oCorrFunc_1D_stream << global_path << "/correlfunct1D" << "_" << temp_particle_name << ".dat";
	ofstream oCorrFunc_1D;
	//oCorrFunc_1D.open(oCorrFunc_1D_stream.str().c_str(), ios::app);
	oCorrFunc_1D.open(oCorrFunc_1D_stream.str().c_str());

	for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
	for(int ipphi = 0; ipphi < n_interp_pphi_pts; ipphi++)
	for(int i = 0; i < qnpts; ++i)
		oCorrFunc_1D << scientific << setprecision(7) << setw(15)
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << q_pts[i] << "   "
			<< CFvals[ipt][ipphi][i][0] << "   " << CFvals[ipt][ipphi][i][1] << "   " << CFvals[ipt][ipphi][i][2] << endl;

	return;
}

void CorrelationFunction::Readin_results(int folderindex)
{
double dummy;
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_ev" << folderindex << no_df_stem << ".dat";
	ifstream inputHBT(filename_stream_HBT.str().c_str());

for(int iKT = 0; iKT < n_localp_T; iKT++)
{
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		inputHBT >> dummy;
		inputHBT >> dummy;
        	inputHBT >> dummy;
		inputHBT >> R2_side[iKT][iKphi];
		inputHBT >> R2_out[iKT][iKphi];
		inputHBT >> R2_outside[iKT][iKphi];
		inputHBT >> R2_long[iKT][iKphi];
		inputHBT >> R2_sidelong[iKT][iKphi];
		inputHBT >> R2_outlong[iKT][iKphi];
	}
}

	inputHBT.close();

	return;
}


void CorrelationFunction::Output_all_dN_dypTdpTdphi(int folderindex)
{
	ostringstream filename_stream_all_dN_dypTdpTdphi;
	filename_stream_all_dN_dypTdpTdphi << global_path << "/all_res_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_all_dN_dypTdpTdphi(filename_stream_all_dN_dypTdpTdphi.str().c_str());
	for(int ii = 0; ii < Nparticle; ii++)
	for(int iphi = 0; iphi < n_interp_pphi_pts; iphi++)
	{
		for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
			output_all_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << spectra[ii][ipt][iphi] << "   ";
		output_all_dN_dypTdpTdphi << endl;
	}
	output_all_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Output_total_target_dN_dypTdpTdphi(int folderindex)
{
	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << global_path << "/total_" << local_name << "_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	for(int iphi = 0; iphi < n_interp_pphi_pts; iphi++)
	{
		for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
			output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << spectra[target_particle_id][ipt][iphi] << "   ";
		output_target_dN_dypTdpTdphi << endl;
	}

	output_target_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Output_chosen_resonances()
{
	ostringstream filename_stream_crf;
	filename_stream_crf << global_path << "/chosen_resonances.dat";
	ofstream output_crf(filename_stream_crf.str().c_str());

	output_crf << particle_monval << endl;
	for (int icr = 0; icr < (int)chosen_resonances.size(); icr++)
		output_crf << all_particles[chosen_resonances[icr]].monval << endl;

	output_crf.close();

	return;
}

/*void CorrelationFunction::Read_in_all_dN_dypTdpTdphi(int folderindex)
{
	for(int wfi = 0; wfi < n_weighting_functions; wfi++)
	{
		ostringstream filename_stream_all_dN_dypTdpTdphi;
		filename_stream_all_dN_dypTdpTdphi << global_path << "/all_res_dN_dypTdpTdphi_mom_"
								<< setfill('0') << setw(2) << wfi << "_ev" << folderindex << no_df_stem << ".dat";
		ifstream input_all_dN_dypTdpTdphi(filename_stream_all_dN_dypTdpTdphi.str().c_str());
	
		int local_filelength = get_filelength(filename_stream_all_dN_dypTdpTdphi.str().c_str());
		int local_filewidth = get_filewidth(filename_stream_all_dN_dypTdpTdphi.str().c_str());
		if (VERBOSE > 0) *global_out_stream_ptr << "Read_in_all_dN_dypTdpTdphi(): nrows = "
							<< local_filelength << " and ncols = " << local_filewidth << endl;
		if ((Nparticle * n_interp_pphi_pts != local_filelength) || (n_interp_pT_pts != local_filewidth))
		{
			cerr << "Read_in_all_dN_dypTdpTdphi(): Mismatch in dimensions in file "
				<< "all_res_dN_dypTdpTdphi_mom_" << setfill('0') << setw(2) << wfi
				<< "_ev" << folderindex << no_df_stem << ".dat!" << endl;
			exit(1);
		}
	
		for(int ii = 0; ii < Nparticle; ii++)
		for(int iphi = 0; iphi < n_interp_pphi_pts; iphi++)
		for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
		{
			input_all_dN_dypTdpTdphi >> dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi];
			if (abs(dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi]) > 1.e-100)
			{
				ln_dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi] = log(abs(dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi]));
				sign_of_dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi] = sgn(dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi]);
			}
			//cout << "Read in (pT, pphi, EdNdp3 ST moms) = " << SPinterp_pT[ipt] << "   " << SPinterp_pphi[iphi] << "   "
			//	<< scientific << setprecision(8) << setw(12) << dN_dypTdpTdphi_moments[ii][wfi][ipt][iphi] << endl;
		}
	
		input_all_dN_dypTdpTdphi.close();
		if (VERBOSE > 0) *global_out_stream_ptr << "Successfully read in " << filename_stream_all_dN_dypTdpTdphi.str().c_str() << endl;
	}


	ostringstream filename_stream_pTpts;
	filename_stream_pTpts << global_path << "/pT_gauss_table.dat";
	ifstream input_pTpts(filename_stream_pTpts.str().c_str());
	ostringstream filename_stream_pphipts;
	filename_stream_pphipts << global_path << "/phi_gauss_table.dat";
	ifstream input_pphipts(filename_stream_pphipts.str().c_str());

	double * dummy_pT_wts = new double [n_interp_pT_pts];
	double * dummy_pphi_wts = new double [n_interp_pphi_pts];

	for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
		input_pTpts >> SPinterp_pT[ipt] >> dummy_pT_wts[ipt];

	for(int ipphi = 0; ipphi < n_interp_pphi_pts; ipphi++)
		input_pphipts >> SPinterp_pphi[ipphi] >> dummy_pphi_wts[ipphi];

	if (VERBOSE > 0) *global_out_stream_ptr << "Read_in_all_dN_dypTdpTdphi(): read in pT and pphi points!" << endl;

	input_pTpts.close();
	input_pphipts.close();

	return;
}*/

//End of file
