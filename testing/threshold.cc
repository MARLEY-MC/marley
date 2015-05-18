// This program reads in the energies of final-state particles from MARLEY
// and rejects those with energy less than 0.5 MeV.
// The output is a .root file containing histograms.


#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TLatex.h"

#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <vector>

int main()
{

  std::ifstream infile("out.txt");
  std::string line; // For holding a line of input
  
  // Histogram to hold measured energies
  TH1D *energy_hist = new TH1D("energy_hist", "Total Measurable Energy", 100, 0., 6.); // MeV
 
  std::vector<double> energy;
  double e, esum;
  bool firstblank = false;
  int count = 0;
  const double threshold = 0.9;

  std::regex empty_line("\\s*");
  std::regex electron("e- kinetic energy = [0-9]+.[0-9]+e[+-][0-9]+");
  std::regex potassium("40K kinetic energy = [0-9]+.[0-9]+e[+-][0-9]+");
  std::regex gamma("gamma energy = [0-9]+.[0-9]+e[+-][0-9]+");
  std::regex electronMass("e- mass = [0-9]+.[0-9]+e[+-][0-9]+");
  std::regex deltaMass("Ground state nuclear mass change = [0-9]+.[0-9]+e[+-][0-9]+");
  
  while(infile)
    {
      e = 0;
      std::getline(infile, line);

      if(!std::regex_match(line, empty_line)) // Separating blocks
	{ 
	  if(std::regex_match(line, gamma))
	    {
	      e = std::stod(line.substr(15,35));
	      std::cout << "Gamma ";
	      firstblank = true;

	      if( e > threshold )
		{
		  std::cout << "above threshold: " << e << std::endl;
		  energy.push_back(e);
		}
	    }
	  else if(std::regex_match(line, electron))
	    {
	      e = std::stod(line.substr(20,40));
	      std::cout << "Electron " << e << std::endl;
	      if( e > threshold )
		{
		  std::cout << "above threshold: " << e << std::endl;
		  energy.push_back(e);
		}
	      count++; // Counting the number of events here, since each event has an e-
	    }
	  
	  else if(std::regex_match(line, electronMass))
	    {
	      e = std::stod(line.substr(10,30));
	      std::cout << "e- mass: " << e << std::endl;
	      energy.push_back(e);
	    }
	  
	  else if(std::regex_match(line, potassium))
	    {
	      e = std::stod(line.substr(21,41));
	      std::cout << "40K " << e << std:: endl;
	      energy.push_back(e);
	    }
	  
	   else if(std::regex_match(line, deltaMass))
	    {
	      e = std::stod(line.substr(35,55));
	      std::cout << "nuclear mass change:  " << e << std::endl;
	      energy.push_back(e);
	    }


	  // Don't record the energy if it's below the threshold
	  // else if (e < 1.0 && e != 0)
	  // std::cout << "below threshold: " << e << std::endl;

	    
	}

      else if ((std::regex_match(line, empty_line) && firstblank))
	  {
	  esum = 0;
	  
	  for(unsigned int i = 0; i < energy.size(); i++)
	    esum += energy[i];

	  energy_hist->Fill(esum);
	  energy.clear();
	  std::cout << "Total energy = " << esum << std::endl << std::endl;

	  firstblank = false;
	}
    }

  std::cout << "Number of entries: " << count << std::endl;
  energy_hist->GetXaxis()->SetTitle("E_{tot} [MeV]");
  
  TFile *rootFile = new TFile("energy_hist.root","RECREATE");
  energy_hist->Write();
  rootFile->Close();  
}
