#include <array>

#include "Particle.hh"

namespace marley {

  /// @brief Simple rotation matrix implementation used to reorient
  /// Particle objects based on the incident neutrino direction
  class RotationMatrix {

    using ThreeVector = std::array<double, 3>;
    using ThreeThreeMatrix = std::array<ThreeVector, 3>;

    public:

      /// @brief Creates a 3&times;3 identity matrix
      RotationMatrix();

      /// @brief Create a 3&times;3 rotation matrix that rotates the 3-vector
      /// from_vec into the 3-vector to_vec
      RotationMatrix(const ThreeVector& from_vec, const ThreeVector& to_vec);

      /// @brief Returns a copy of the 3-vector v normalized to have unit
      /// magnitude
      static ThreeVector normalize(const ThreeVector& v);

      /// @brief Create a rotated copy of the 3-vector v
      ThreeVector rotate_copy(const ThreeVector& v);

      /// @brief Rotate a 3-vector v in place
      void rotate_inplace(ThreeVector& v);

      /// @brief Rotate the 3-momentum of a marley::Particle in place
      void rotate_particle_inplace(marley::Particle& p);

    protected:

      /// @brief 3&times;3 rotation matrix
      ThreeThreeMatrix matrix_;
  };

}
