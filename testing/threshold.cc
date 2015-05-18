// This program reads in the energies of final-state particles from MARLEY
// and rejects those with energy less than 0.5 MeV.


#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"

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
  TH1D *energy_hist = new TH1D("energy_hist", "Total Energy", 100, 3., 7.); // MeV
 
  std::vector<double> energy;
  double e, esum;
  bool firstblank = false;

  std::regex empty_line("\\s*");
  std::regex electron("e- kinetic energy = [0-9]+.[0-9]+e[+-][0-9]+");
  std::regex potassium("40K kinetic energy = [0-9]+.[0-9]+e[+-][0-9]+");
  std::regex gamma("gamma energy = [0-9]+.[0-9]+e[+-][0-9]+");
  
  while(infile)
    {
      e = 0;
      std::getline(infile, line);
      //std::cout << line << std::endl;

      if(!std::regex_match(line, empty_line)) // Separating blocks
	{ 
	  if(std::regex_match(line, gamma))
	    {
	      e = std::stod(line.substr(15,35));
	      std::cout << "Gamma "; // << e << std::endl;
	      firstblank = true;
	    }
	  else if(std::regex_match(line, electron))
	    {
	      e = std::stod(line.substr(20,40));
	      std::cout << "Electron "; // << e << std::endl;
	    }
	  else if(std::regex_match(line, potassium))
	    {
	      e = std::stod(line.substr(21,41));
	      std::cout << "40K "; // << e << std::endl;
	    }

	  // else
	  //std::cout << line << std::endl;
	      
	  if( e > 0.5 ) // Threshold of 0.5 MeV
	    {
	      std::cout << "above threshold: " << e << std::endl;
	      energy.push_back(e);
	    }

	  // Don't record the energy if it's below the threshold
	   else if (e < 0.5 && e != 0)
	     std::cout << "below threshold: " << e << std::endl;

	    
	}
      
      else if (std::regex_match(line, empty_line) && firstblank)
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

  TFile *rootFile = new TFile("energy_hist.root","RECREATE");
  energy_hist->Write();
  rootFile->Close();  
}
