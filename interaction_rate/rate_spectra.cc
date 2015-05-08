//  This class can be used to generate an approximate supernova 
//  neutrino interaction spectrum in a liquid argon detector. The
//  supernova distance is assumed to be 10 kpc.  The user can select
//  the CC reaction (electron neutrino or electron antineutrino) as 
//  shown by the examples below:
//
//  ./rate_spectra antineutrino
//
//  or
//
//  ./rate_spectra neutrino
//
//
//  C. Grant.
//------------------------------------------------------------ 

#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <string>

float *sn_neutrino_flux(std::string name, float energy[], float flux[]){

  // Find neutrino flux file
  if(name.compare("neutrino") == 0){
    std::ifstream nu_e_file("../sn_neutrino_flux/sn_electron_neutrino_flux.txt");
    if(nu_e_file.is_open()){
      // Format is:  "Neutrino Energy [MeV]" "Number neutrinos"
      for(int i=0; i<100; i++){
	nu_e_file >> energy[i] >> flux[i];
      }
    }
    else std::cout << "Unable to open electron neutrino flux file." << std::endl;
    if(nu_e_file.is_open()) nu_e_file.close();
  }


  // Find antineutrino flux file
  if(name.compare("antineutrino") == 0){
    std::ifstream anti_nu_e_file("../sn_neutrino_flux/sn_electron_antineutrino_flux.txt");
    if(anti_nu_e_file.is_open()){
      // Format is:  "Neutrino Energy [MeV]" "Number neutrinos"
      for(int i=0; i<100; i++){
	anti_nu_e_file >> energy[i] >> flux[i];
      }
    }
    else std::cout << "Unable to open electron anti neutrino flux file" << std::endl;
    if(anti_nu_e_file.is_open()) anti_nu_e_file.close();
  }

  return flux;
}

float *total_cross_section(std::string name, float xsection[]){
  float energy;

  // Find neutrino flux file
  if(name.compare("neutrino") == 0){
    std::ifstream nu_e_file("../total_cross_section/kolbe_electron_neutrino_absorption.txt");
    if(nu_e_file.is_open()){
      // Format is:  "Neutrino Energy [MeV]" "Number neutrinos"
      for(int i=0; i<100; i++){
	nu_e_file >> energy >> xsection[i];
      }
    }
    else std::cout << "Unable to open electron neutrino xsection file." << std::endl;
    if(nu_e_file.is_open()) nu_e_file.close();
  }


  // Find antineutrino flux file
  if(name.compare("antineutrino") == 0){
    std::ifstream anti_nu_e_file("../total_cross_section/kolbe_electron_antineutrino_absorption.txt");
    if(anti_nu_e_file.is_open()){
      // Format is:  "Neutrino Energy [MeV]" "Number neutrinos"
      for(int i=0; i<100; i++){
	anti_nu_e_file >> energy >> xsection[i];
      }
    }
    else std::cout << "Unable to open electron anti neutrino xsection file" << std::endl;
    if(anti_nu_e_file.is_open()) anti_nu_e_file.close();
  }

  return xsection;
}


int main(int argc, char *argv[]){

  std::string neutrino_type = argv[1];
  
  std::cout << " " << std::endl;
  std::cout << "Obtaining flux file for " << neutrino_type << " type..." << std::endl;
  std::cout << " " << std::endl;

  float energy[100];
  float flux[100];

  sn_neutrino_flux(neutrino_type,energy,flux);

  std::cout << " " << std::endl;
  std::cout << "Obtaining cross-section file for " << neutrino_type << " type..." << std::endl;
  std::cout << " " << std::endl;

  float xsection[100];

  total_cross_section(neutrino_type,xsection);

  std::cout << " " << std::endl;
  std::cout << "Producing interaction spectrum for " << neutrino_type << " type..." << std::endl;
  std::cout << " " << std::endl;

  std::string rate_file_name = "rate_spectrum_electron_";
  rate_file_name += neutrino_type;
  rate_file_name += ".txt";

  float rate[100] = {0};
  float total_rate = 0;

  std::ofstream rate_spectrum(rate_file_name.c_str());
  if(rate_spectrum.is_open()){
    for(int i=0; i<100; i++){
      rate[i] = flux[i]*xsection[i];
      total_rate += rate[i];
    }
    
    for(int i=0; i<100; i++){
      rate_spectrum << energy[i] << "\t" << rate[i]/total_rate << "\n";
    }
  }
  else std::cout << "Unable to open rate file" << std::endl;
  if(rate_spectrum.is_open()) rate_spectrum.close();

  std::cout << " " << std::endl;
  std::cout << "Finished!" << std::endl;
  std::cout << " " << std::endl;



}
