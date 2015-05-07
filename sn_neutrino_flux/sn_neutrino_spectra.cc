//  This class can be used to generate an approximate supernova 
//  neutrino distribution for use as input.  The spectrum comes
//  from a type II supernova releasing 3x10^53 ergs of energy.
//  It is based on the Fermi-Dirac distribution described in 
//  K. Langanke, et al. PRL 76 (1996) 2629 and arXiv:hep-ph/0307222v1.
//  The total cross-section as a function of energy is folded 
//  with the Fermi-Dirac distribution to obtain the number of events.
//
//
//  C. Grant.
//------------------------------------------------------------ 

#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>

double fermi_dirac_distribution(float C, bool e_flavor, bool anti, float nu_energy){
  float eta = 0;
  float T = 0;
  float N_nu;

  if(e_flavor && !anti){
    T = 3.5; // temperature in MeV
    N_nu = 2.8;  // total number of electron neutrinos expected (x10^57)
  }
  else if(e_flavor && anti){
    T = 5.0; // temperature in MeV
    N_nu = 1.9;  // total number of electron anti-neutrinos expected (x10^57)
  }
  else if(!e_flavor && (!anti || anti)){
    T = 8.0; // temperature in MeV
    N_nu = 5.0; // total number of mu+tau neutrinos + anti-neutrinos expected (x10^57)
  }
  return (C/std::pow(T,3))*(std::pow(nu_energy,2)/(1+std::exp(nu_energy/(T-eta))))*N_nu;
}


int main(){
  
  std::cout << "Calculating flux files..." << std::endl;

  // First perform the electron flavor calculation for C = 0.55
  std::ofstream nu_e_file("sn_electron_neutrino_flux.txt");
  if(nu_e_file.is_open()){
   // Format is:  "Neutrino Energy [MeV]" "Number neutrinos"
    for(int i=0; i<100; i++){
      nu_e_file << i << "\t" << fermi_dirac_distribution(0.55,true,false,i) << "\n";
    }
  }
  else std::cout << "Unable to open electron neutrino flux file." << std::endl;
  nu_e_file.close();


  // Next perform the anti electron flavor calculation for C = 0.55
  std::ofstream anti_nu_e_file("sn_electron_antineutrino_flux.txt");
  if(anti_nu_e_file.is_open()){
    // Format is:  "Neutrino Energy [MeV]" "Number neutrinos"
    for(int i=0; i<100; i++){
      anti_nu_e_file << i << "\t" << fermi_dirac_distribution(0.55,true,true,i) << "\n";
    }
  }
  else std::cout << "Unable to open electron anti neutrino flux file" << std::endl;
  anti_nu_e_file.close();


  // Finally perform all other flavor calculation for C = 0.55
  std::ofstream other_nu_file("sn_other_neutrino_flux.txt");
  if(other_nu_file.is_open()){
    // Format is:  "Neutrino Energy [MeV]" "Number neutrinos"
    for(int i=0; i<100; i++){
      other_nu_file << i << "\t" << fermi_dirac_distribution(0.55,false,true,i) << "\n";
    }
  }
  else std::cout << "Unable to open other neutrino flux file" << std::endl;
  other_nu_file.close();

  std::cout << "Finished!" << std::endl;

}
