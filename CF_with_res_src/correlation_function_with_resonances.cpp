#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<sys/time.h>
#include<algorithm>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include "../src/Stopwatch.h"
#include "../src/parameters.h"
#include "../src/readindata.h"
#include "../src/CFWR.h"
#include "../src/generate_processing_record.h"
#include "../src/plumberglib.h"
#include "../src/sorter.h"
#include "correlation_function_with_resonances.h"

using namespace std;

int main(int argc, char *argv[])
{
   Stopwatch sw;
   Stopwatch sw_total;
   sw_total.tic();
   sw.tic();

   bool generatedcorrfuncs = false;
   string currentworkingdirectory = get_selfpath();
   int folderindex = get_folder_index(currentworkingdirectory);
   initialize_PRfile(currentworkingdirectory);

   ostringstream filename_stream;
   filename_stream << currentworkingdirectory << "/Processing_record.txt";
   ofstream output(filename_stream.str().c_str(), ios::app);

   output << "/**********Processing output**********/" << endl;
   output << "entering folder: " << currentworkingdirectory << endl;

   //load freeze out and particle information
   int FO_length = 0;
   int particle_idx = 1;  //for pion+

   ostringstream decdatfile;
   output << "Loading the decoupling data...." << endl;
   decdatfile << currentworkingdirectory << "/decdat2.dat";
   output << decdatfile.str() << endl;
   FO_length=get_filelength(decdatfile.str().c_str());
   output << "Total number of freeze out fluid cell: " <<  FO_length << endl;

   //read the data arrays for the decoupling information
   FO_surf* FOsurf_ptr = new FO_surf[FO_length];
   read_decdat(FO_length, FOsurf_ptr, currentworkingdirectory, false);
   
   //read the positions of the freeze out surface
   read_surfdat(FO_length, FOsurf_ptr, currentworkingdirectory);
   
   //read the chemical potential on the freeze out surface
   int N_stableparticle;
   ifstream particletable("/home/plumberg.1/HBTPlumberg/EOS/EOS_particletable.dat");
   particletable >> N_stableparticle;
   double** particle_mu = new double* [N_stableparticle];
   for(int i=0; i<N_stableparticle; i++)
      particle_mu[i] = new double [FO_length];
   for(int i=0; i<N_stableparticle; i++)
      for(int j=0; j<FO_length; j++)
         particle_mu[i][j] = 0.0;
   if(N_stableparticle >0)
   {
      //if(hydropara_ptr->IEOS==7)       //for s95p_PCE
         read_decdat_mu(FO_length, N_stableparticle, particle_mu, currentworkingdirectory);
   }

   //read particle resonance decay table
   particle_info *particle = new particle_info [Maxparticle];
   int Nparticle=read_resonance(particle);
   output <<"read in total " << Nparticle << " particles!" << endl;
   output << "Calculating "<< particle[particle_idx].name << endl;
   if(N_stableparticle >0)
   {
      output << " EOS is partially chemical equilibrium " << endl;
      calculate_particle_mu(7, Nparticle, FOsurf_ptr, FO_length, particle, particle_mu);
   }
   else
   {
	output << " EOS is chemical equilibrium. " << endl;
	for(int j=0; j<Nparticle; j++)
	{
		particle[j].mu = 0.0e0;
		for(int i=0; i<FO_length; i++)
			FOsurf_ptr[i].particle_mu[j] = 0.0e0;
	}
   }
	//calculate (semi-analytic approximation of) pure thermal spectra for all particle species
	calculate_thermal_particle_yield(Nparticle, particle, FOsurf_ptr[0].Tdec);
	//use this to estimate resonance-decay contributions from each particles species to final state particle, here, pion(+),
	//as well as set effective branching ratios
	compute_total_contribution_percentages(particle_idx, Nparticle, particle);
	//sort all particles by importance of their percentage contributions, then compute resonance SVs for only contributions up to some threshold
   sw.toc();
   output << "read in data finished!" << endl;
   output << "Used " << sw.takeTime() << " sec." << endl;

   //HBT calculations begin ...
   double localy = 0.0e0;
   sw.tic();
     if(fabs(localy) > 1e-16)
     {
        output << "Case of y != 0 not yet supported.  Exiting..." << endl;
        return 0;
     }

	//double threshold = 0.35;	//include only enough of the most important resonances to account for fixed fraction of total resonance-decay pion(+)s
	double threshold = atof(argv[1]);
				//threshold = 1.0 means include all resonance-decay pion(+)s,
				//threshold = 0.0 means include none of them
	double net_fraction_resonance_contribution = 0.0;
	output << "Working with threshold = " << threshold << endl;
	vector<int> chosen_resonance_indices;
	if (threshold > 1.0 + 1.e-10)
	{
		int single_chosen_resonance = 183;	//index in all_particles array
		chosen_resonance_indices.push_back(single_chosen_resonance);
		output << "\t --> Looking only at contributions from " << particle[single_chosen_resonance].name << " and descendants" << endl;
		//get_all_descendants(&chosen_resonance_indices, particle, Nparticle, output);
		sort_by_mass(&chosen_resonance_indices, particle, Nparticle, output);
	}
	else
	{
		get_important_resonances(particle_idx, &chosen_resonance_indices, particle, Nparticle, threshold, net_fraction_resonance_contribution, output);
		get_all_descendants(&chosen_resonance_indices, particle, Nparticle, output);
		sort_by_mass(&chosen_resonance_indices, particle, Nparticle, output);
		for (int ii = 0; ii < (int)chosen_resonance_indices.size(); ii++)
			output << ii << "   " << chosen_resonance_indices[ii] << "   " << particle[chosen_resonance_indices[ii]].name << endl;
	}

   CorrelationFunction correlation_function(&particle[particle_idx], particle, Nparticle, FOsurf_ptr, chosen_resonance_indices, particle_idx, output);
   correlation_function.read_in_all_dN_dypTdpTdphi = false;
   correlation_function.output_all_dN_dypTdpTdphi = !(correlation_function.read_in_all_dN_dypTdpTdphi);
   correlation_function.currentfolderindex = folderindex;
   correlation_function.Set_path(currentworkingdirectory);
   correlation_function.Set_use_delta_f(true);
   correlation_function.Set_ofstream(output);

   bool rescale_truncated_resonance_contributions = true;
   if (rescale_truncated_resonance_contributions)
	correlation_function.fraction_of_resonances = net_fraction_resonance_contribution;
   else
	correlation_function.fraction_of_resonances = 1.0;	//i.e., no rescaling

   correlation_function.Update_sourcefunction(&particle[particle_idx], FO_length, particle_idx);

   output << "Calculating HBT radii via Gaussian fit method..." << endl;
   correlation_function.Compute_correlation_function(FOsurf_ptr);
   correlation_function.Get_GF_HBTradii(FOsurf_ptr, folderindex);	//does outputting of results too
   
   correlation_function.Output_total_target_dN_dypTdpTdphi(folderindex);
   correlation_function.Output_chosen_resonances();
   correlation_function.Output_results(folderindex);
   output << "Finished calculating HBT radii via Gaussian fit method" << endl;



   sw.toc();
   output << "Finished in " << sw.takeTime() << " sec." << endl;
   sw_total.toc();
   output << "Program totally finished in " << sw_total.takeTime() << " sec." << endl;

   output << "/**********End of processing output**********/" << endl;

   output.close();

   checkforfiles_PRfile(currentworkingdirectory, folderindex, generatedcorrfuncs);

   finalize_PRfile(currentworkingdirectory);

   return 0;
}
