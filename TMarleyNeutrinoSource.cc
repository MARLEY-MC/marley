#include "TMarleyGenerator.hh"
#include "TMarleyNeutrinoSource.hh"

TMarleyNeutrinoSource::TMarleyNeutrinoSource(double Emin, double Emax,
  NeutrinoType nu_type)
{
  type = nu_type;
  E_min = Emin;
  E_max = Emax;

  if(type == NeutrinoType::ElectronNeutrino) {
    temperature = 3.5;
    // total number of electron neutrinos expected (x10^57)
    tot_num_nu = 2.8;
  }
  else if(type == NeutrinoType::ElectronAntineutrino) {
    temperature = 5.0;
    // total number of electron anti-neutrinos expected (x10^57)
    tot_num_nu = 1.9;
  }
  else { // muon or tauon flavor
    temperature = 8.0;
    // total number of mu+tau neutrinos + anti-neutrinos expected (x10^57)
    tot_num_nu = 5.0;
  }

  // Create a forwarding call wrapper for Fermi-Dirac distribution function that
  // takes a single argument. This will be used for rejection sampling of
  // neutrino energies from this source.
  fd_dist = std::bind(&fermi_dirac_distribution, temperature, tot_num_nu,
    std::placeholders::_1 /*nu_energy*/);
}

double TMarleyNeutrinoSource::fermi_dirac_distribution(double T, double N_nu,
  double nu_energy)
{
  return (C / std::pow(T, 3)) * (std::pow(nu_energy, 2)
    / (1 + std::exp(nu_energy / (T - eta)))) * N_nu;
}

double TMarleyNeutrinoSource::sample_neutrino_energy(TMarleyGenerator& gen) {
  return gen.rejection_sample(fd_dist, E_min, E_max);
}
