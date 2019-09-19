#pragma once
#include <string>

#include "marley/Event.hh"
#include "marley/MassTable.hh"

namespace marley {

  class Generator;

  /// @brief Abstract base class that represents a two-two scattering reaction
  /// @details This class models a reaction of the form a + b &rarr; c + d.
  /// The projectile (particle a) is taken to have lab-frame kinetic energy KEa
  /// and to be traveling toward the target along the positive z direction.
  /// The target (particle b) is taken to be at rest in the lab frame.
  /// @todo Add check that projectile kinetic energy KEa >= 0. for all
  /// relevant member functions of Reaction and NuclearReaction
  class Reaction {

    public:

      virtual ~Reaction() = default;

      // Enumerated type describing the kind of scattering process
      // represented by a Reaction
      enum class ProcessType { CC, NC };

      /// @brief Compute the reaction's total cross section (MeV<sup> -2</sup>)
      /// @param pdg_a Projectile's PDG code
      /// @param KEa Lab-frame kinetic energy of the incident projectile
      /// @return %Reaction total cross section (MeV<sup> -2</sup>)
      /// @note Functions that override total_xs() should always return 0.
      /// if pdg_a != pdg_a_.
      virtual double total_xs(int pdg_a, double KEa) = 0;

      /// @brief Differential cross section
      /// @f$d\sigma/d\cos\theta_{c}^{\mathrm{CM}}@f$
      /// (MeV<sup> -2</sup>)
      /// @param pdg_a PDG code for the projectile
      /// @param KEa Lab-frame kinetic energy (MeV) of the projectile
      /// @param cos_theta_c_cm CM frame scattering cosine of the ejectile
      /// @note This function should return 0. if pdg_a != pdg_a_.
      virtual double diff_xs(int pdg_a, double KEa,
        double cos_theta_c_cm) = 0;

      /// @brief Create an event object for this reaction
      /// @param pdg_a PDG code for the incident projectile
      /// @param KEa Lab-frame kinetic energy of the projectile
      /// @param gen Reference to the Generator to use for random sampling
      virtual marley::Event create_event(int pdg_a, double KEa,
        marley::Generator& gen) = 0;

      /// @brief Get a string that contains the formula for this reaction
      inline std::string get_description() { return description_; }

      /// @brief Get the process type for this reaction
      inline ProcessType process_type() { return process_type_; }

    protected:

      int pdg_a_; ///< PDG code for the projectile
      int pdg_b_; ///< PDG code for the target
      int pdg_c_; ///< PDG code for the ejectile
      int pdg_d_; ///< PDG code for the residue

      double ma_; ///< Projectile mass (MeV)
      double mb_; ///< Target mass (MeV)
      double mc_; ///< Ejectile mass (MeV)
      double md_; ///< Residue mass (MeV)

      /// @brief String that contains a formula describing the reaction
      std::string description_;

      /// @brief Type of scattering process (CC, NC) represented by
      /// this reaction
      ProcessType process_type_;

      /// @brief Helper function that handles CM frame kinematics for the
      /// reaction
      /// @param KEa Lab-frame kinetic energy (MeV) of the projectile
      /// @param[out] s <a
      /// href="https://en.wikipedia.org/wiki/Mandelstam_variables"> Mandelstam
      /// s</a> (MeV<sup>2</sup>)
      /// @param[out] Ec_cm Ejectile total energy (MeV) in the CM frame
      /// @param[out] pc_cm Ejectile 3-momentum magnitude (MeV) in the CM frame
      /// @param[out] Ed_cm Residue total energy (MeV) in the CM frame
      void two_two_scatter(double KEa, double& s, double& Ec_cm,
        double& pc_cm, double& Ed_cm);

      /// @brief Helper function that makes an event object.
      /// @details This function should be called by
      /// marley::Reaction::create_event() after CM frame scattering angles
      /// have been sampled for the ejectile. For reactions where the residue
      /// may be left in an excited state, the excitation energy should be
      /// recorded by supplying it as the final argument.
      /// @param KEa Lab-frame kinetic energy (MeV) of the projectile
      /// @param pc_cm Ejectile 3-momentum magnitude (MeV) in the CM frame
      /// @param cos_theta_c_cm Cosine of ejectile's CM frame polar angle
      /// @param phi_c_cm Ejectile's CM frame azimuthal angle (radians)
      /// @param Ec_cm Ejectile total energy (MeV) in the CM frame
      /// @param Ed_cm Residue total energy (MeV) in the CM frame
      /// @param E_level Residue excitation energy (MeV)
      virtual marley::Event make_event_object(double KEa,
        double pc_cm, double cos_theta_c_cm, double phi_c_cm,
        double Ec_cm, double Ed_cm, double E_level = 0.);
  };

}
