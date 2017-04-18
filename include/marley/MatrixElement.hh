#pragma once

namespace marley {

  /// @brief A nuclear matrix element that represents a transition
  /// caused by a neutrino-nucleus reaction
  class MatrixElement {

    public:

      /// @brief Enumerated type that represents the possible kinds of nuclear
      /// transitions recognized by MARLEY
      enum TransitionType { FERMI = 0, GAMOW_TELLER = 1 };

      /// @param level_energy Excitation energy (MeV) of the final-state nuclear
      /// level accessed by the matrix element
      /// @param strength Numerical value (dimensionless) of the matrix element
      /// @param type Type of nuclear transition (e.g., Fermi, Gamow-Teller)
      /// represented by the matrix element
      /// @param final_level Pointer to the Level object that represents the
      /// final nuclear level
      inline MatrixElement(double level_energy, double strength, TransitionType
        type, marley::Level* final_level = nullptr)
        : level_energy_(level_energy), strength_(strength), type_(type),
        final_level_(final_level) {}

      /// @brief Get the excitation energy (MeV) of the final-state nuclear
      /// level accessed by the matrix element
      inline double level_energy() const { return level_energy_; }

      /// @brief Get the numerical value (dimensionless) of the matrix element
      inline double strength() const { return strength_; }

      /// @brief Get the kind of nuclear transition (e.g., Fermi, Gamow-Teller)
      /// represented by the matrix element
      inline TransitionType type() const { return type_; }

      /// @brief Get a pointer to the final-state nuclear Level accessed by the
      /// matrix element, or nullptr if it is a transition to the unbound
      /// continuum
      inline const marley::Level* level() const { return final_level_; }

      /// @brief Get a pointer to the final-state nuclear Level accessed by the
      /// matrix element, or nullptr if it is a transition to the unbound
      /// continuum
      inline marley::Level* level() { return final_level_; }

      /// @brief Set the excitation energy (MeV) of the final-state nuclear
      /// level accessed by the matrix element
      inline void set_level_energy(double energy) { level_energy_ = energy; }

      /// @brief Set the numerical value (dimensionless) of the matrix element
      inline void set_strength(double strength) { strength_ = strength; }

      /// @brief Set the kind of nuclear transition (e.g., Fermi, Gamow-Teller)
      /// represented by the matrix element
      inline void set_type(TransitionType type) { type_ = type; }

      /// @brief Set the pointer to the final nuclear Level object accessed
      /// by the matrix element
      inline void set_level(marley::Level* lev) { final_level_ = lev; }

    protected:

      /// @brief Energy (MeV) of the final-state nuclear level accessed by this
      /// matrix element
      double level_energy_;

      /// @brief Numerical value of the matrix element (dimensionless)
      double strength_;

      /// @brief The kind of transition represented by this matrix element
      /// (Fermi, Gamow-Teller, etc.)
      TransitionType type_;

      /// @brief Pointer to the final Level object for a transition to a
      /// discrete nuclear level, or nullptr for an unbound state
      marley::Level* final_level_;
  };
}
