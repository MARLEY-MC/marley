#pragma once
#include <functional>
#include <regex>
#include <string>
#include <vector>

#include "DecayScheme.hh"
#include "Event.hh"
#include "Level.hh"
#include "MassTable.hh"
#include "Reaction.hh"
#include "StructureDatabase.hh"

namespace marley {

  class Generator;

  class NuclearReaction : public marley::Reaction {
    public:
      NuclearReaction(std::string filename, marley::StructureDatabase& db);
      double fermi_function(double beta_c);
      double fermi_approx(int Z, double E, bool electron);
      double max_level_energy(double Ea);
      double get_threshold_energy();
      // Total reaction cross section (in MeV^(-2)), including all final nuclear
      // levels, for an incident projectile with lab-frame total energy Ea
      virtual double total_xs(int pdg_a, double Ea);
      void set_decay_scheme(marley::DecayScheme* scheme);
      marley::Event create_event(int particle_id_a, double Ea,
        marley::Generator& gen);
      double sample_cos_theta_c_cm(/*double matrix_el,*/ int m_type,
        double beta_c_cm, marley::Generator& gen);

      inline size_t get_num_levels() const {
        return residue_level_energies.size();
      }

      inline double get_level_energy(size_t index, bool& bound,
        double& strength)
      {
        strength = residue_level_strengths[index];

        // Get the excitation energy for the current level
        if (residue_level_pointers[index] != nullptr) {
          bound = true;
          return residue_level_pointers[index]->get_energy();
        }
        else {
          bound = false;
          // Level is unbound, so just use the energy given in the reaction dataset
          // rather than trying to find its ENSDF version
          return residue_level_energies[index];
        }
      }

      inline virtual marley::Event make_event_object(double Ea, double pc_cm,
        double cos_theta_c_cm, double phi_c_cm, double Ec_cm, double Ed_cm,
        double E_level = 0.)
      {
        marley::Event event = marley::Reaction::make_event_object(Ea, pc_cm,
          cos_theta_c_cm, phi_c_cm, Ec_cm, Ed_cm, E_level);
        // Assume that the target is a neutral atom (q_b = 0)
        event.target().set_charge(0);
        // Assign the correct charge to the residue
        event.residue().set_charge(q_d);
        // Return the completed event object
        return event;
      }

    private:
      // Residue ground state mass
      double md_gs;
      // Initial and final values of the atomic and mass numbers
      int Zi, Ai, Zf, Af;
      // Net charge of particle d (in units of the proton charge) following this reaction
      int q_d;

      // Total cross section for a given final nuclear level energy, in units
      // convenient for sampling
      virtual double total_xs(double E_level, double Ea, double matrix_element);

      // Lab-frame total energy of the projectile at threshold for this reaction
      // (all final-state particles at rest in the CM frame)
      double Ea_threshold;

      std::vector<double> residue_level_energies; // Energy values from reaction dataset
      std::vector<double> residue_level_strengths; // B(F) and B(GT) values from reaction dataset
      // TODO: come up with a more general way of representing what kind of matrix element
      // is given in a reaction dataset entry. Right now 0 represents B(F) and 1 represents B(GT),
      // but when you start to include other kinds of matrix elements (e.g., forbidden transitions)
      // this will need to change.
      std::vector<double> residue_level_strength_ids; // Matrix element type identifiers from reaction dataset
      std::vector<marley::Level*> residue_level_pointers; // Pointers to the corresponding ENSDF decay scheme levels
  };

}
