#pragma once
#include "marley_utils.hh"

class TMarleyGenerator;

class TMarleyNeutrinoSource {
  public:
    enum class NeutrinoType { ElectronNeutrino, MuonNeutrino, TauonNeutrino,
      ElectronAntineutrino, MuonAntineutrino, TauonAntineutrino };

    // TODO: improve or remove this default constructor
    // Use the default neutrino source settings
    inline TMarleyNeutrinoSource() {
      type = NeutrinoType::ElectronNeutrino;
      E_min = 4.36;
      E_max = 50;
      temperature = 3.5;
      tot_num_nu = 2.8;
      fd_dist = std::bind(&fermi_dirac_distribution, temperature, tot_num_nu,
        std::placeholders::_1 /*nu_energy*/);
    }

    TMarleyNeutrinoSource(double Emin, double Emax, NeutrinoType nu_type
      = NeutrinoType::ElectronNeutrino);

    // Sample a neutrino energy using the Fermi-Dirac distribution
    double sample_neutrino_energy(TMarleyGenerator& gen);

    inline NeutrinoType get_neutrino_type() {
      return type;
    }

  private:
    // Toy model of supernova neutrino spectrum
    static double fermi_dirac_distribution(double T, double N_nu,
      double nu_energy);
    // Normalization constant for Fermi-Dirac distribution
    static constexpr double C = 0.55;
    // Temperature offset for Fermi-Dirac distribution
    static constexpr double eta = 0;
    // Temperature for Fermi-Dirac distribution (in MeV)
    double temperature;
    // Total number of neutrinos expected (x10^57)
    double tot_num_nu;
    // Neutrino type produced by this source object
    NeutrinoType type;
    // Minimum and maximum neutrino energies produced by this source
    double E_min, E_max;
    // Function wrapper used for sampling from the Fermi-Dirac distribution
    std::function<double(double)> fd_dist;
};
