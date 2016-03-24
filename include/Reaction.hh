#pragma once
#include <string>

#include "Event.hh"
#include "MassTable.hh"

namespace marley {

  class Generator;
  
  // Abstract base class that represents a 2-2 scattering reaction
  // between  a projectile (particle a) with lab-frame total energy Ea,
  // and a target (particle b) that is taken to be at rest in the lab frame.
  // Particles a and b scatter into particles c (the ejectile) and d (the
  // residue).
  class Reaction {
    public:
      // Total reaction cross section (in MeV^(-2)) for an incident
      // projectile with lab-frame total energy Ea
      virtual double total_xs(int particle_id_a, double Ea) = 0;
      // Creates an event object for this reaction using the generator gen
      virtual marley::Event create_event(int particle_id_a, double Ea,
        marley::Generator& gen) = 0;
  
      inline std::string get_description() { return description; }
  
    protected:
      // Particle ID numbers (PDG convention)
      int pid_a, pid_b, pid_c, pid_d;
      // Particle masses and squared masses (pre-computed for speed)
      double ma, mb, mc, md;
      double ma2, mb2, mc2, md2;
  
      // String describing the reaction
      std::string description;
  
      // Handles kinematics for the reaction, loading variables with the value of
      // Mandelstam s, the CM frame ejectile total energy, the CM frame ejectile
      // 3-momentum magnitude, and the CM frame residue total energy
      void two_two_scatter(double Ea, double& s, double& Ec_cm,
        double& pc_cm, double& Ed_cm);
  
      // Helper function that makes an event object. This should be called
      // by create_event after CM frame scattering angles have been sampled
      // for the ejectile. For reactions where the residue may be left in an
      // excited state, the excitation energy may be recorded by supplying
      // it via the final argument.
      virtual marley::Event make_event_object(double Ea,
        double pc_cm, double cos_theta_c_cm, double phi_c_cm,
        double Ec_cm, double Ed_cm, double E_level = 0.);
  };

}
