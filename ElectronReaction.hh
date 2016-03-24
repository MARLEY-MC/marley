#pragma once
#include <functional>
#include <string>

#include "Event.hh"
#include "MassTable.hh"
#include "Reaction.hh"

namespace marley {

  class Generator;
  
  // Represents a reaction in which an incident (anti)neutrino of any flavor
  // elastically scatters off of an electron originally bound to an atom with
  // atomic number Zatom. The cross section computed for this reaction is
  // computed per atom (scaled up from the single-electron cross section by the
  // atomic number Zatom, which is the same as the number of electron targets
  // bound to the atom).
  class ElectronReaction : public marley::Reaction {
    public:
      ElectronReaction(size_t Z);
  
      // Total reaction cross section (in MeV^(-2)) for an incident
      // projectile with lab-frame total energy Ea
      virtual double total_xs(int particle_id_a, double Ea);
  
      // Creates an event object for this reaction using the generator gen
      virtual marley::Event create_event(int particle_id_a, double Ea,
        marley::Generator& gen);
  
    private:
      // Atomic number of the atom used as the target by this reaction
      size_t Zatom;
      // Coupling constants to use for cross section calculations (updated
      // each time based on the projectile's particle ID)
      double g1, g2;
      // Helper variable for g2 (pre-computed for speed)
      double g2_squared_over_three;
  
      // Helper function for cross section calculations. Returns true if the
      // projectile's particle ID was recognized, or false if it was not. Loads
      // g1 and g2 with the appropriate coupling constants.
      bool determine_coupling_constants(int particle_id_a);
  
      // CM frame differential cross section in units convenient for sampling.
      // The input variables are Mandelstam s and the ejectile's CM frame
      // scattering cosine.
      double diff_xs_cm_for_sampling(double s, double cos_theta_c_cm);
  
      // Sample a CM frame scattering cosine for the ejectile. The first argument
      // is Mandelstam s.  This function assumes that the coupling constants g1
      // and g2 have already been updated, so make sure to use it only after
      // calling determine_coupling_constants.
      double sample_cos_theta_c_cm(double s, marley::Generator& gen);
  };

}
