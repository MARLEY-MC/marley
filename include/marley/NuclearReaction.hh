#pragma once
#include <functional>
#include <regex>
#include <string>
#include <vector>

#include "marley/DecayScheme.hh"
#include "marley/Event.hh"
#include "marley/Level.hh"
#include "marley/MassTable.hh"
#include "marley/Reaction.hh"
#include "marley/StructureDatabase.hh"

namespace marley {

  class Generator;

  /// @brief A neutrino-nucleus reaction
  class NuclearReaction : public marley::Reaction {

    public:

      /// @param filename Name (with path, if needed) of the matrix element
      /// data file
      /// @param db Reference to the StructureDatabase to use for sampling
      /// excited nuclear levels
      NuclearReaction(std::string filename, marley::StructureDatabase& db);

      virtual marley::Event create_event(int particle_id_a,
        double KEa, marley::Generator& gen) override;

      /// @brief Compute the
      /// <a href="https://en.wikipedia.org/wiki/Beta_decay#Fermi_function">
      /// Fermi function</a>
      /// @param beta_c <a
      /// href="http://scienceworld.wolfram.com/physics/RelativisticBeta.html">
      /// Dimensionless speed</a> of the ejectile
      double fermi_function(double beta_c) const;

      /// @brief Get the maximum possible excitation energy (MeV) of the
      /// final-state residue that is kinematically allowed
      /// @param KEa Projectile lab-frame kinetic energy (MeV)
      double max_level_energy(double KEa) const;

      /// @brief Sets the DecayScheme object to use for sampling excited levels
      /// in the residue
      void set_decay_scheme(marley::DecayScheme* scheme);

      /// @brief Get the minimum lab-frame kinetic energy (MeV) of the
      /// projectile that allows this reaction to proceed via a transition to
      /// the residue's ground state
      double threshold_kinetic_energy() const;

      /// @brief Total reaction cross section (MeV<sup> -2</sup>), including
      /// all kinematically-allowed final nuclear levels
      /// @param pdg_a PDG code for the projectile
      /// @param KEa Lab-frame kinetic energy (MeV) of the projectile
      /// @return %Reaction total cross section (MeV<sup> -2</sup>)
      /// @note This function returns 0. if pdg_a != pdg_a_.
      virtual double total_xs(int pdg_a, double KEa) override;

    private:

      virtual marley::Event make_event_object(double KEa,
        double pc_cm, double cos_theta_c_cm, double phi_c_cm, double Ec_cm,
        double Ed_cm, double E_level = 0.) override;

      /// @brief Samples a polar angle cosine for the ejectile using
      /// the relevant portion of the reaction nuclear matrix element
      /// @param m_type Integer representing the type of transition
      /// (0 = Fermi, 1 = Gamow-Teller)
      /// @param beta_c_cm <a
      /// href="http://scienceworld.wolfram.com/physics/RelativisticBeta.html">
      /// Dimensionless speed</a> of the ejectile in the CM frame
      /// @param gen Reference to the Generator to use for random sampling
      double sample_cos_theta_c_cm(/*double matrix_el,*/ int m_type,
        double beta_c_cm, marley::Generator& gen) const;

      /// @brief Total cross section for a given final nuclear level energy, in
      /// units convenient for sampling
      /// @param E_level Residue excitation energy (MeV)
      /// @param KEa Lab-frame projectile kinetic energy (MeV)
      /// @param matrix_element Nuclear matrix element for a transition
      /// to the level of interest
      virtual double total_xs(double E_level, double KEa,
        double matrix_element) const;

      double md_gs_; ///< Ground state mass (MeV) of the residue

      int Zi_; ///< Target atomic number
      int Ai_; ///< Target mass number
      int Zf_; ///< Residue atomic number
      int Af_; ///< Residue mass number

      /// @brief Net charge of the residue (in units of the proton charge)
      /// following this reaction
      int q_d_;

      /// @brief Lab-frame kinetic energy of the projectile at threshold for
      /// this reaction (i.e., the residue is produced in its ground state, and
      /// all final-state particles are at rest in the CM frame)
      double KEa_threshold_;

      /// @todo Refactor NuclearReaction class to use a vector of objects
      /// rather than multiple vectors to hold the level energies, strengths,
      /// etc.

      /// @brief Excitation energy values (MeV) for each level from the
      /// reaction data file
      std::vector<double> residue_level_energies_;

      /// @brief B(F) and B(GT) values from the reaction data file
      std::vector<double> residue_level_strengths_;

      /// @brief Matrix element type identifiers from the reaction data file
      /// @todo Come up with a more general way of representing what kind of
      /// matrix element is given in a reaction dataset entry. Right now 0
      /// represents B(F) and 1 represents B(GT), but when you start to include
      /// other kinds of matrix elements (e.g., forbidden transitions) this
      /// will need to change. One possible solution is an enum or enum class
      /// that can be extended as needed.
      std::vector<double> residue_level_strength_ids_;

      /// @brief Pointers to the Level objects representing final states
      /// of the residue.
      /// @details Particle-unbound levels have nullptr entries and de-excite
      /// using the HauserFeshbachDecay model.
      std::vector<marley::Level*> residue_level_pointers_;
  };

}
