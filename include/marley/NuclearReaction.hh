#pragma once
#include <functional>
#include <regex>
#include <string>
#include <vector>

#include "marley/DecayScheme.hh"
#include "marley/Event.hh"
#include "marley/Level.hh"
#include "marley/MassTable.hh"
#include "marley/MatrixElement.hh"
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
      NuclearReaction(const std::string& filename,
        marley::StructureDatabase& db);

      /// Produces a two-two scattering Event that proceeds via this reaction
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

      /// @brief Differential cross section
      /// @f$d\sigma/d\cos\theta_{c}^{\mathrm{CM}}@f$
      /// (MeV<sup> -2</sup>) summed over all kinematically-allowed final
      /// nuclear levels
      /// @param pdg_a PDG code for the projectile
      /// @param KEa Lab-frame kinetic energy (MeV) of the projectile
      /// @param cos_theta_c_cm CM frame scattering cosine of the ejectile
      /// @note This function returns 0. if pdg_a != pdg_a_.
      virtual double summed_diff_xs(int pdg_a, double KEa,
        double cos_theta_c_cm);

      /// @brief Total cross section (MeV<sup> -2</sup>) for a given final
      /// nuclear level
      /// @param E_level Residue excitation energy (MeV)
      /// @param KEa Lab-frame projectile kinetic energy (MeV)
      /// @param matrix_element Nuclear matrix element for a transition
      /// to the level of interest
      /// @param tr_type Enum value indicating what kind of nuclear transition
      /// (e.g., Fermi or Gamow-Teller) is represented by the matrix element
      /// @param[out] beta_c_cm After this function is called, this variable
      /// will contain the (dimensionless) speed of the ejectile as measured
      /// in the center of momentum frame
      /// @param check_max_E_level Whether to check the reaction threshold
      /// via a call to max_level_energy(). If this argument is set to false,
      /// then the check will be skipped.
      virtual double total_xs(double E_level, double KEa,
        double matrix_element, marley::MatrixElement::TransitionType tr_type,
        double& beta_c_cm, bool check_max_E_level = true) const;

      /// Computes an approximate correction factor to account for
      /// effects of the Coulomb potential when calculating cross sections
      /// @param beta_rel_cd The relative speed of the final particles c and d
      /// (dimensionless)
      double coulomb_correction_factor(double beta_rel_cd) const;

      /// Computes a Coulomb correction factor according to the effective
      /// momentum approximation. See J. Engel, Phys. Rev. C 57, 2004 (1998)
      /// @param beta_rel_cd The relative speed of the final particles c and d
      /// (dimensionless)
      /// @param ok Flag that is set to false if subtracting the Coulomb
      /// potential pulls the event below threshold
      double ema_factor(double beta_rel_cd, bool& ok) const;

      /// Computes the weak nuclear charge @f$ Q_W = N - (1
      /// - 4\sin^2\theta_W)Z @f$ for the target nucleus
      /// @details In the expression above, @f$ N @f$ (@f$ Z @f$) is the
      /// neutron (proton) number of the target nucleus and
      /// @f$ \theta_W @f$ is the weak mixing angle.
      double weak_nuclear_charge() const;

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
      double sample_cos_theta_c_cm(const marley::MatrixElement& matrix_el,
        double beta_c_cm, marley::Generator& gen) const;

      /// Helper function for total_xs and summed_diff_xs()
      /// @param pdg_a PDG code for the projectile
      /// @param KEa Lab-frame projectile kinetic energy (MeV)
      /// @param cos_theta_c_cm CM frame scattering cosine for the ejectile
      /// @param differential Whether the total cross section (false)
      /// or the differential cross section (true,
      /// @f$d\sigma/d\cos\theta_{c}^{\mathrm{CM}}@f$) should
      /// be summed over all kinematically accessible nuclear levels
      /// by this function
      /// @return The requested cross section (MeV<sup> -2</sup>)
      double summed_xs_helper(int pdg_a, double KEa, double cos_theta_c_cm,
        bool differential);

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

      /// @brief Matrix elements representing all of the possible nuclear
      /// transitions that may be caused by this reaction
      std::vector<marley::MatrixElement> matrix_elements_;
  };

}
