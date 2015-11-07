#pragma once
#include <functional>
#include <regex>
#include <string>
#include <vector>

#include "TMarleyDecayScheme.hh"
#include "TMarleyEvent.hh"
#include "TMarleyLevel.hh"
#include "TMarleyMassTable.hh"
#include "TMarleyStructureDatabase.hh"

class TMarleyGenerator;

class TMarleyReaction {
  public:
    TMarleyReaction(std::string filename, TMarleyStructureDatabase& db);
    double fermi_function(int Z, int A, double E, bool electron);
    double fermi_approx(int Z, double E, bool electron);
    double max_level_energy(double Ea);
    double get_threshold_energy();
    // Total reaction cross section (in MeV^(-2)), including all final nuclear levels
    double total_xs_cm(double Ea);
    // Total cross section for a given final nuclear level energy, in
    // units convenient for sampling
    double total_xs_cm(double E_level, double Ea, double matrix_element);
    void set_decay_scheme(TMarleyDecayScheme* scheme);
    TMarleyEvent create_event(double Ea, TMarleyGenerator& gen);
    double sample_cos_theta_c_cm(/*double matrix_el,*/ int m_type,
      double beta_c_cm, TMarleyGenerator& gen);
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

  private:
    double ma; // projectile mass
    double mb; // target mass
    double mc; // ejectile mass
    double md_gs; // residue mass
    double ma2, mb2, mc2; // Squared masses
    static constexpr double GF = 1.16637e-11; // Fermi coupling constant (MeV^(-2)) 
    static constexpr double Vud = 0.97427; // abs(V_ud) (from CKM matrix)
    int Zi, Ai, Zf, Af; // Initial and final values of the atomic and mass numbers
    int pid_a, pid_b, pid_c, pid_d; // Particle IDs for all 4 particles

    // Lab-frame total energy of the projectile at threshold
    // for this reaction (all final-state particles at rest
    // in the CM frame)
    double Ea_threshold;

    std::vector<double> residue_level_energies; // Energy values from reaction dataset
    std::vector<double> residue_level_strengths; // B(F) and B(GT) values from reaction dataset
    // TODO: come up with a more general way of representing what kind of matrix element
    // is given in a reaction dataset entry. Right now 0 represents B(F) and 1 represents B(GT),
    // but when you start to include other kinds of matrix elements (e.g., forbidden transitions)
    // this will need to change.
    std::vector<double> residue_level_strength_ids; // Matrix element type identifiers from reaction dataset
    std::vector<TMarleyLevel*> residue_level_pointers; // Pointers to the corresponding ENSDF decay scheme levels
};
