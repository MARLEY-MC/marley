/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/LICENSE or
// visit http://opensource.org/licenses/GPL-3.0

#include <array>
#include <cmath>
#include <string>

#include "marley/Error.hh"
#include "marley/RotationMatrix.hh"

using ThreeVector = std::array<double, 3>;

// Anonymous namespace for helper functions
namespace {

  // Loads the 3-vector dest with the cross product v1 x v2
  void cross_product(ThreeVector& dest,
    const ThreeVector& v1, const ThreeVector& v2)
  {
    dest[0] = v1[1] * v2[2] - v1[2] * v2[1];
    dest[1] = v1[2] * v2[0] - v1[0] * v2[2];
    dest[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }

  // Returns the dot product v1 . v2
  double dot_product(const ThreeVector& v1,
    const ThreeVector& v2)
  {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  }

  // Loads the 3-vector dest with the difference v1 - v2
  void subtract(ThreeVector& dest,
    const ThreeVector& v1, const ThreeVector& v2)
  {
    dest[0] = v1[0] - v2[0];
    dest[1] = v1[1] - v2[1];
    dest[2] = v1[2] - v2[2];
  }

}

marley::RotationMatrix::RotationMatrix()
  : matrix_{{ {{ 1., 0., 0.}}, {{ 0., 1., 0.}}, {{ 0., 0., 1.}} }}
{}

// Returns a copy of the 3-vector v normalized to have unit magnitude
ThreeVector marley::RotationMatrix::normalize(const ThreeVector& v)
{
  static ThreeVector nv({0., 0., 0.});
  double norm_factor = std::sqrt(std::pow(v[0], 2) + std::pow(v[1], 2)
    + std::pow(v[2], 2));
  if (norm_factor <= 0.) throw marley::Error(std::string("Invalid vector")
    + " magnitude encountered in marley::RotationMatrix::normalize()");
  else norm_factor = 1. / norm_factor;
  nv[0] = norm_factor * v[0];
  nv[1] = norm_factor * v[1];
  nv[2] = norm_factor * v[2];
  return nv;
}

// Returns a rotated copy of the 3-vector v
ThreeVector marley::RotationMatrix::rotate_copy(const ThreeVector& v)
{
  ThreeVector rv = {0., 0., 0.};
  for (unsigned i = 0; i < 3; ++i) rv[i] = dot_product(matrix_[i], v);
  return rv;
}

// Rotates a 3-vector v in place
void marley::RotationMatrix::rotate_inplace(ThreeVector& v)
{
  ThreeVector rv = {0., 0., 0.};
  for (unsigned i = 0; i < 3; ++i) rv[i] = dot_product(matrix_[i], v);
  v = rv;
}

// Rotates the 3-momentum of a marley::Particle in place
void marley::RotationMatrix::rotate_particle_inplace(marley::Particle& p)
{
  ThreeVector rv = {0., 0., 0.};
  ThreeVector three_momentum = { p.px(), p.py(), p.pz() };

  for (unsigned i = 0; i < 3; ++i)
    rv[i] = dot_product(matrix_[i], three_momentum);

  p.set_px(rv[0]);
  p.set_py(rv[1]);
  p.set_pz(rv[2]);
}

/// @details <p>This function is a C++11 version of an original rotation matrix
/// program by M&ouml;ller &amp; Hughes (see
/// <a href="http://tinyurl.com/hperc7d">this</a> GitHub page for details)</p>
/// <p>The vectors from_vec and to_vec do not need to be normalized, but both
/// should be nonzero. If either vector is a null vector, then a marley::Error
/// will be thrown.</p>
marley::RotationMatrix::RotationMatrix(const ThreeVector& from_vec,
  const ThreeVector& to_vec)
{
  static constexpr ThreeVector null_three_vector = { 0., 0., 0. };

  if (from_vec == null_three_vector)
    throw marley::Error(std::string("Null from vector")
      + " passed to constructor of marley::RotationMatrix");
  else if (to_vec == null_three_vector)
    throw marley::Error(std::string("Null to vector")
      + " passed to constructor of marley::RotationMatrix");

  // The from and to vectors must be normalized to use this method, so make any
  // necessary adjustments now.
  ThreeVector from = normalize(from_vec);
  ThreeVector to = normalize(to_vec);

  double e = dot_product(from, to);
  double f = std::abs(e);

  static constexpr double EPSILON = 0.000001;
  if (f > 1.0 - EPSILON) { // "from" and "to" vectors are almost parallel

    // Temporary storage vectors
    ThreeVector v;
    ThreeVector u;

    // Find the standard unit vector x most nearly orthogonal to "from"
    ThreeVector x;
    x[0] = std::abs(from[0]);
    x[1] = std::abs(from[1]);
    x[2] = std::abs(from[2]);

    if (x[0] < x[1])
    {
      if (x[0] < x[2])
      {
        x[0] = 1.0; x[1] = x[2] = 0.0;
      }
      else
      {
        x[2] = 1.0; x[0] = x[1] = 0.0;
      }
    }
    else
    {
      if (x[1] < x[2])
      {
        x[1] = 1.0; x[0] = x[2] = 0.0;
      }
      else
      {
        x[2] = 1.0; x[0] = x[1] = 0.0;
      }
    }

    // u = x - from
    subtract(u, x, from);

    // v = x - to
    subtract(v, x, to);

    // coefficients for later use
    double c1 = 2.0 / dot_product(u, u);
    double c2 = 2.0 / dot_product(v, v);
    double c3 = c1 * c2  * dot_product(u, v);

    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        matrix_[i][j] =  - c1 * u[i] * u[j] - c2 * v[i] * v[j]
          + c3 * v[i] * u[j];
      }
      matrix_[i][i] += 1.0;
    }
  }
  else  // the most common case, unless "from" = "to", or "from" = -"to"
  {
    ThreeVector v;
    cross_product(v, from, to); // v = from x to

    // hand-optimized version (9 mults less than original)
    // optimization by Gottfried Chen
    double hvx, hvz, hvxy, hvxz, hvyz;
    double h = 1.0 / (1.0 + e);
    hvx = h * v[0];
    hvz = h * v[2];
    hvxy = hvx * v[1];
    hvxz = hvx * v[2];
    hvyz = hvz * v[1];

    matrix_[0][0] = e + hvx * v[0];
    matrix_[0][1] = hvxy - v[2];
    matrix_[0][2] = hvxz + v[1];

    matrix_[1][0] = hvxy + v[2];
    matrix_[1][1] = e + h * v[1] * v[1];
    matrix_[1][2] = hvyz - v[0];

    matrix_[2][0] = hvxz - v[1];
    matrix_[2][1] = hvyz + v[0];
    matrix_[2][2] = e + hvz * v[2];
  }
}
