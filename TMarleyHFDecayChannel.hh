#pragma once

#include "TMarleyFragment.hh"
#include "TMarleyLevel.hh"
#include "TMarleyParity.hh"

class TMarleyHFDecayChannel {
  public:
    inline const TMarleyFragment& get_fragment() const {
      return fragment;
    }

    // Determine final values for several quantities of
    // interest for the residual nuclide after the decay occurs
    virtual void get_post_decay_parameters(double& Exf, int& two_Jf,
      TMarleyParity& Pi_f) = 0;

  protected:
    const TMarleyFragment& fragment;
};

class TMarleyHFDiscreteDecayChannel : public TMarleyHFDecayChannel {
  public:

    inline TMarleyHFDiscreteDecayChannel(const TMarleyFragment& frag,
      TMarleyLevel* f_lev): fragment(frag), final_level(f_lev)
    {
    }

    inline TMarleyLevel* get_final_level() {
      return final_level;
    }

    inline void get_post_decay_parameters(double& Exf, int& two_Jf,
      TMarleyParity& Pi_f)
    {
      Exf = final_level->get_energy();
      two_Jf = final_level->get_two_J();
      Pi_f = final_level->get_parity();
    }

  private:
    TMarleyLevel* final_level; // Pointer to final level in a nuclear decay scheme,
                               // nullptr if in the unbound continuum
    //int l;      // Orbital angular momentum of the outgoing fragment
};

class TMarleyHFContinuumDecayChannel : public TMarleyHFDecayChannel {
  public:
    inline TMarleyHFContinuumDecayChannel(const TMarleyFragment& frag,
      TMarleyGenerator& gener, std::function<double(double)>& cpw_func,
      const TMarleySphericalOpticalModel& om,
      double Emin, double Emax, double width):
      fragment(frag), gen(gener), om(optical_model)
    {
      E_c_min = Emin;
      Exf_max = Emax;
      total_width = width;

      cpw = [&total_width, cpw_func] (double energy) -> double
        { return cpw_func(energy) / total_width; };
    }

    void get_post_decay_parameters(double& Exf, int& two_Jf,
      TMarleyParity& Pi_f)
    {
      Exf = gen.rejection_sample(cpw, double E_c_min, double Exf_max);

      get_continuum_jpi_widths(Exf);
      std::discrete_distribution<size_t>::param_type params(j_pi_widths.begin(),
        j_pi_widths.end());
      size_t jpi_index = gen.discrete_sample(j_pi_dist, params);

      two_Jf = two_Js.at(jpi_index);
      Pi_f = Pis.at(jpi_index);
    }

    inline void get_fragment_continuum_jpi_widths(double Exf) {

      two_Js.clear();
      Pis.clear();
      j_pi_widths.clear();

      TMarleyParity Pa = fragment.get_parity();
      int two_s = fragment.get_two_s();
      int fragment_pid = fragment.get_pid();
    
      double Ea = (Mconst - Exf*(2*Mfgs_ion + Exf)) / (2 * Mi);
      int Zf = om.get_Z();
      int Af = om.get_A();
      // Final nuclear state parity
      TMarleyParity Pf;
      // The orbital parity starts as (-1)^0 = 1. Rather than applying parity
      // conservation each time, just find the final state parity Pf for l = 0.
      // Then we can safely flip Pf without further thought for each new l value in
      // the loop.
      if (Pi == Pa) Pf = 1;
      else Pf = -1;
      // For each new iteration, increment l and flip the final-state parity
      for (int l = 0; l <= l_max; ++l, !Pf) {
        int two_l = 2*l;
        for (int two_j = std::abs(two_l - two_s);
          two_j <= two_l + two_s; two_j += 2)
        {
          for (int twoJf = std::abs(twoJi - two_j);
            twoJf <= twoJi + two_j; twoJf += 2)
          {
            double width = om.transmission_coefficient(Ea, fragment_pid, two_j, l,
              two_s, DEFAULT_NUMEROV_STEP_SIZE)
              * TMarleyBackshiftedFermiGasModel::level_density(Zf, Af, Exf, twoJf);
    
            two_Js.push_back(twoJf);
            Pis.push_back(Pif)
            j_pi_widths.push_back(width);
          }
        }
      }
    }

  private:
    double E_c_min; // Minimum energy for the continuum
    double Exf_max; // Maximum accessible final excitation energy
    double total_width; // Total width to the continuum, used for normalizing the PDF when
                        // sampling the final nuclear excitation energy
    double Mconst;
    double Mfgs_ion;
    double Mi;
    TMarleyGenerator& gen;
    const TMarleySphericalOpticalModel& om;
    // Normalized continuum partial decay width as a function of energy
    std::function<double(double)> cpw;
    //int l;      // Orbital angular momentum of the outgoing fragment
    std::vector<double> j_pi_widths;
    std::vector<int> two_Js;
    std::vector<TMarleyParity> Pis;
    std::discrete_distribution<size_t> j_pi_dist;
};
