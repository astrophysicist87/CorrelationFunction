#ifndef GENERATE_PROCESSING_RECORD_H
#define GENERATE_PROCESSING_RECORD_H

#include <string>
#include <sstream>
#include <time.h>

using namespace std;

#include "plumberglib.h"
#include "parameters.h"

void initialize_PRfile(string currentworkingdirectory, string PRfilename = "Processing_record.txt")
{
	ostringstream filename_stream;
	filename_stream << currentworkingdirectory << "/" << PRfilename;
	ofstream output(filename_stream.str().c_str(), ios::app);

	string lambdaflagstring = lambdaflag ? truestring : falsestring;

	output << "/***************************************************/" << endl;
	output << "/****************" << PRfilename << "**************/" << endl;
	output << "/***************************************************/" << endl;


	output << "General initializations:" << endl;
	output << "   - Spatial rapidity information:" << endl;
	output << "      --> eta_s_npts: " << eta_s_npts << endl;
	output << "      --> eta_s_i: " << eta_s_i << endl;
	output << "      --> eta_s_f: " << eta_s_f << endl;

	output << "   - Relative momentum information:" << endl;
	output << "      --> qnpts: " << qnpts << endl;
	output << "      --> delta_q: " << delta_q << endl;
	output << "      --> init_q: " << init_q << endl;

	output << "   - Correlation function information:" << endl;
	output << "      --> corrfuncdim: " << corrfuncdim << endl;
	output << "      --> lambdaflag: " << lambdaflagstring << endl;

	output << "   - Pair momentum information:" << endl;
	output << "      --> n_localp_T: " << n_localp_T << endl;
	output << "      --> localp_T_min: " << localp_T_min << endl;
	output << "      --> localp_T_max: " << localp_T_max << endl;
	output << "      --> n_localp_phi: " << n_localp_phi << endl;
	output << "      --> localp_phi_min: " << localp_phi_min << endl;
	output << "      --> localp_phi_max: " << localp_phi_max << endl;

	output << "   - HBT Fourier information:" << endl;
	output << "      --> n_order: " << n_order << endl;

	output << "   - Miscellaneous information:" << endl;
	output << "      --> CWD: " << currentworkingdirectory << endl;

	output << "/***************************************************/" << endl << endl;

	output.close();

	return;
}

void checkforfiles_PRfile(string currentworkingdirectory, int folderindex, bool corrfuncsgenerated, string PRfilename = "Processing_record.txt")
{
	ostringstream filename_stream;
	filename_stream << currentworkingdirectory << "/" << PRfilename;
	ofstream output(filename_stream.str().c_str(), ios::app);

	stringstream HBTSVfilename, HBTcfsSVfilename, planepsifilename;
	stringstream HBTGFfilename, HBTcfsGFfilename;
	HBTSVfilename << currentworkingdirectory << "/HBTradii_ev" << folderindex << ".dat" << endl;
	HBTcfsSVfilename << currentworkingdirectory << "/HBTradii_cfs_ev" << folderindex << ".dat" << endl;
	HBTGFfilename << currentworkingdirectory << "/HBTradii_GF_ev" << folderindex << ".dat" << endl;
	HBTcfsGFfilename << currentworkingdirectory << "/HBTradii_GF_cfs_ev" << folderindex << ".dat" << endl;
	planepsifilename << currentworkingdirectory << "/plane_psi_ev" << folderindex << ".dat" << endl;

	string HBTSVexists = fexists(HBTSVfilename.str().c_str()) ? truestring : falsestring;
	string HBTcfsSVexists = fexists(HBTcfsSVfilename.str().c_str()) ? truestring : falsestring;
	string HBTGFexists = fexists(HBTGFfilename.str().c_str()) ? truestring : falsestring;
	string HBTcfsGFexists = fexists(HBTcfsGFfilename.str().c_str()) ? truestring : falsestring;
	string planepsiexists = fexists(planepsifilename.str().c_str()) ? truestring : falsestring;
	string corrfuncsexist = corrfuncsgenerated ? truestring : falsestring;

	output << "/***************************************************/" << endl
		<< "HBT output files (source variances method):" << endl
		<< "      --> Generated " << HBTSVfilename.str() << ": " << HBTSVexists << endl
		<< "      --> Generated " << HBTcfsSVfilename.str() << ": " << HBTcfsSVexists << endl
		<< "      --> Generated " << planepsifilename.str() << ": " << planepsiexists << endl;

	output << "HBT output files (Gaussian fit method):" << endl
		<< "      --> Generated " << HBTGFfilename.str() << ": " << HBTGFexists << endl
		<< "      --> Generated " << HBTcfsGFfilename.str() << ": " << HBTcfsGFexists << endl;

	output << "Correlation function output files:" << endl
		<< "      --> Generated correlation functions in corrfuncs_ev" << folderindex << ".zip: " << corrfuncsexist << endl;

	output.close();

	return;
}

void finalize_PRfile(string currentworkingdirectory, string PRfilename = "Processing_record.txt")
{
	ostringstream filename_stream;
	filename_stream << currentworkingdirectory << "/" << PRfilename;
	ofstream output(filename_stream.str().c_str(), ios::app);

	stringstream out;

	time_t now = time(0);
	tm *ltm = localtime(&now);

	char* dt = ctime(&now);

	out << setfill('0') << setw(2) << 1 + ltm->tm_mon << "/" << setw(2) << ltm->tm_mday << "/" << setw(4) << 1900 + ltm->tm_year;

	string date = out.str();

	output << "Timestamp: " << dt << endl;

	output << "/***************************************************/" << endl;
	output << "/*********End of processing for " << date << "**********/" << endl;
	output << "/***************************************************/" << endl;

	output.close();

	return;
}

#endif
