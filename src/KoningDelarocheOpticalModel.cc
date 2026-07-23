/// @file
/// @copyright Copyright (C) 2016-2023 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

#include "marley/marley_utils.hh"
#include "marley/HauserFeshbachDecay.hh"
#include "marley/Logger.hh"
#include "marley/KoningDelarocheOpticalModel.hh"

#include "marley/coulomb_wavefunctions.hh"

std::complex<double>
marley::KoningDelarocheOpticalModel::optical_model_potential(double r,
  double fragment_KE_lab, int fragment_pdg, int two_j, int l, int two_s,
  int target_charge)
{
  const double target_mass = target_mass_for_charge( target_charge );

  // Kinematic-state construction needs the total CM energy, so retrieve the
  // fragment mass directly before constructing the per-request state.
  /// @todo TODO: find a better way of doing this!
  const auto& mt = marley::MassTable::Instance();
  double m_fragment = mt.get_particle_mass( fragment_pdg );

  double KE_tot_CM = std::max(0., marley_utils::real_sqrt(
    std::pow(target_mass + m_fragment, 2)
    + 2.*target_mass*fragment_KE_lab) - m_fragment - target_mass);

  auto state = make_working_state( KE_tot_CM, fragment_pdg, target_mass );
  calculate_om_parameters(fragment_pdg, two_j, l, two_s, state);
  return omp(r, state);
}

// Finish an optical model potential calculation by taking the r
// dependence into account. Don't add in the Coulomb potential.
std::complex<double>
marley::KoningDelarocheOpticalModel::omp_minus_Vc(double r,
  const WorkingState& state) const
{
  double f_v = f(r, state.Rv, state.av);
  double dfdr_d = dfdr(r, state.Rd, state.ad);

  double temp_Vv = state.Vv * f_v;
  double temp_Wv = state.Wv * f_v;
  double temp_Wd = -4 * state.Wd * state.ad * dfdr_d;

  double temp_Vso = 0;
  double temp_Wso = 0;

  if (state.spin_orbit_eigenvalue != 0) {

    double factor_so = lambda_piplus2 * dfdr(r, state.Rso, state.aso)
      * state.spin_orbit_eigenvalue / r;

    temp_Vso = state.Vso * factor_so;
    temp_Wso = state.Wso * factor_so;
  }

  return std::complex<double>(-temp_Vv + temp_Vso,
    -temp_Wv - temp_Wd + temp_Wso);
}

// Compute all of the pieces of the optical model that depend on the fragment's
// kinetic energy in the lab frame but not on its distance from the origin r.
// Store them in the request-local working state.
void marley::KoningDelarocheOpticalModel::calculate_om_parameters(
  int fragment_pdg, int two_j, int l, int two_s, WorkingState& state) const
{
  // Fragment atomic, mass, and neutron numbers
  state.z = marley_utils::get_particle_Z(fragment_pdg);
  int a = marley_utils::get_particle_A(fragment_pdg);
  int n = a - state.z;

  // Abbreviate the variable name here for simplicity
  const double E = state.fragment_KE_lab;

  // Eigenvalue of the spin-orbit operator
  // 2*(l.s) = j*(j + 1)  - l*(l + 1) -  s*(s + 1)
  // = 0.25*((2j - 2s)*(2j + 2s + 2)) - l*(l+1)
  // (to keep the units right we take hbar = 1).
  bool spin_zero = two_s == 0;
  if (spin_zero) state.spin_orbit_eigenvalue = 0;
  else state.spin_orbit_eigenvalue = 0.25*((two_j - two_s)
    * (two_j + two_s + 2)) - l*(l + 1);


  // Geometrical parameters
  state.Rv = 0;
  state.av = 0;
  state.Rd = 0;
  state.ad = 0;
  state.Rso = 0;
  state.aso = 0;

  // Terms in the spherical optical model potential
  state.Vv = 0;
  state.Wv = 0;
  state.Wd = 0;
  state.Vso = 0;
  state.Wso = 0;

  // Energy to use when computing folded potentials
  double E_eff = E / a;

  if (n > 0) {
    double Ediff_n = E_eff - Efn;
    double Ediff_n2 = std::pow(Ediff_n, 2);
    double Ediff_n3 = std::pow(Ediff_n, 3);

    state.Vv += n * v1n * (1 - v2n*Ediff_n + v3n*Ediff_n2 - v4n*Ediff_n3);
    state.Wv += n * w1n * Ediff_n2 / (Ediff_n2 + std::pow(w2n, 2));
    state.Wd += n * d1n * Ediff_n2 * std::exp(-d2n * Ediff_n)
      / (Ediff_n2 + std::pow(d3n, 2));

    state.Rv += n * Rvn;
    state.av += n * avn;
    state.Rd += n * Rdn;
    state.ad += n * adn;
    state.Rso += n * Rso_n;
    state.aso += n * aso_n;

    if (!spin_zero) {
      double Ediff_so_n = E - Efn;
      double Ediff_so_n2 = std::pow(Ediff_so_n, 2);
      state.Vso += vso1n * std::exp(-vso2n * Ediff_so_n);
      state.Wso += wso1n * Ediff_so_n2
        / (Ediff_so_n2 + std::pow(wso2n, 2));
    }
  }

  if (state.z > 0) {
    double Ediff_p = E_eff - Efp;
    double Ediff_p2 = std::pow(Ediff_p, 2);
    double Ediff_p3 = std::pow(Ediff_p, 3);

    state.Vv += state.z * v1p
      * (1 - v2p*Ediff_p + v3p*Ediff_p2 - v4p*Ediff_p3
      + Vcbar_p*(v2p - 2*v3p*Ediff_p + 3*v4p*Ediff_p2));
    state.Wv += state.z * w1p * Ediff_p2
      / (Ediff_p2 + std::pow(w2p, 2));
    state.Wd += state.z * d1p * Ediff_p2 * std::exp(-d2p * Ediff_p)
      / (Ediff_p2 + std::pow(d3p, 2));

    state.Rv += state.z * Rvp;
    state.av += state.z * avp;
    state.Rd += state.z * Rdp;
    state.ad += state.z * adp;
    state.Rso += state.z * Rso_p;
    state.aso += state.z * aso_p;

    if (!spin_zero) {
      double Ediff_so_p = E - Efp;
      double Ediff_so_p2 = std::pow(Ediff_so_p, 2);
      state.Vso += vso1p * std::exp(-vso2p * Ediff_so_p);
      state.Wso += wso1p * Ediff_so_p2
        / (Ediff_so_p2 + std::pow(wso2p, 2));
    }
  }

  if (a > 1) {
    state.Rv /= a;
    state.av /= a;
    state.Rd /= a;
    state.ad /= a;
    state.Rso /= a;
    state.aso /= a;

    // Apply folding factor for composite particle spin-orbit potentials
    if (!spin_zero) {
      bool z_odd = state.z % 2;
      bool n_odd = n % 2;
      // This factor stays zero for even-even nuclides (which should all be
      // spin-zero anyway)
      double factor = 0.;
      if (z_odd && n_odd) factor = 2.0; // odd-odd
      else if (z_odd != n_odd) factor = 1.0; // even-odd
      factor /= 2*a;

      state.Vso *= factor;
      state.Wso *= factor;
    }
  }
}

marley::KoningDelarocheOpticalModel::KoningDelarocheOpticalModel(int Z,
  int A, double step_size) : marley::OpticalModel(Z, A), step_size_(step_size)
{
  int N = A_ - Z_; // Neutron number

  double A_to_the_one_third = std::pow(A_, 1.0/3.0);

  // Initialize the spherical optical model parameters (see TALYS 1.6 manual)

  // Neutrons
  v1n = 59.30 - 21.0*(N - Z_)/A_ - 0.024*A_; // MeV
  v2n = 0.007228 - 1.48e-6*A_; // MeV^(-1)
  v3n = 1.994e-5 - 2.0e-8*A_; // MeV^(-2)
  v4n = 7e-9; // MeV^(-3)
  w1n = 12.195 + 0.0167*A_; // MeV
  w2n = 73.55 + 0.0795*A_; // MeV
  d1n = 16.0 - 16.0*(N - Z_)/A_; // MeV
  d2n = 0.0180 + 0.003802/(1 + std::exp((A_ - 156.)/8.0)); // MeV^(-1)
  d3n = 11.5; // MeV
  vso1n = 5.922 + 0.0030*A_; // MeV
  vso2n = 0.0040; // MeV^(-1)
  wso1n = -3.1; // MeV
  wso2n = 160.; // MeV
  Efn = -11.2814 + 0.02646*A_; // MeV
  Rvn = 1.3039*A_to_the_one_third - 0.4054; // fm
  avn = 0.6778 - 1.487e-4*A_; // fm
  Rdn = 1.3424*A_to_the_one_third
    - 0.01585*std::pow(A_to_the_one_third, 2); // fm
  adn = 0.5446 - 1.656e-4*A_; // fm
  Rso_n = 1.1854*A_to_the_one_third - 0.647; // fm
  aso_n = 0.59; // fm

  // Protons
  v1p = 59.30 + 21.0*(N - Z_)/A_ - 0.024*A_; // MeV
  v2p = 0.007067 + 4.23e-6*A_; // MeV^(-1)
  v3p = 1.729e-5 + 1.136e-8*A_; // MeV^(-2)
  v4p = v4n; // MeV^(-3)
  w1p = 14.667 + 0.009629*A_; // MeV
  w2p = w2n; // MeV
  d1p = 16.0 + 16.0*(N - Z_)/A_; // MeV
  d2p = d2n; // MeV^(-1)
  d3p = d3n; // MeV
  vso1p = vso1n; // MeV
  vso2p = vso2n; // MeV^(-1)
  wso1p = wso1n; // MeV
  wso2p = wso2n; // MeV
  Efp = -8.4075 + 0.01378*A_; // MeV
  Rvp = Rvn;
  avp = avn;
  Rdp = Rdn;
  adp = 0.5187 + 5.205e-4*A_; // fm
  Rso_p = Rso_n;
  aso_p = aso_n;
  Rc = 1.198*A_to_the_one_third + 0.697/A_to_the_one_third
    + 12.994*std::pow(A_to_the_one_third, -4); // fm
  Vcbar_p = 1.73 * Z_ / Rc; // MeV

}

double marley::KoningDelarocheOpticalModel::total_cross_section(
  double fragment_KE_lab, int fragment_pdg, int two_s, size_t l_max,
  int target_charge)
{
  const double target_mass = target_mass_for_charge( target_charge );

  // Kinematic-state construction needs the total CM energy, so retrieve the
  // fragment mass directly before constructing the per-request state.
  /// @todo TODO: find a better way of doing this!
  const auto& mt = marley::MassTable::Instance();
  double m_fragment = mt.get_particle_mass( fragment_pdg );

  double KE_tot_CM = std::max(0., marley_utils::real_sqrt(
    std::pow(target_mass + m_fragment, 2)
    + 2.*target_mass*fragment_KE_lab) - m_fragment - target_mass);

  const auto state = make_working_state( KE_tot_CM, fragment_pdg,
    target_mass );

  double sum = 0.;
  for (size_t l = 0; l <= l_max; ++l) {
    int two_l = 2*l;
    for (int two_j = std::abs(two_l - two_s);
      two_j <= two_l + two_s; two_j += 2)
    {
      std::complex<double> S = s_matrix_element(fragment_pdg, two_j, l,
        two_s, state);
      sum += (two_j + 1) * (1 - S.real());
    }
  }

  // Compute the cross section in natural units (MeV^(-2))
  double xs = marley_utils::two_pi * sum
    / ((two_s + 1) * state.CM_frame_momentum_squared);
  return xs;
}

double marley::KoningDelarocheOpticalModel::transmission_coefficient(
  double total_KE_CM, int fragment_pdg, int two_j, int l, int two_s,
  int target_charge)
{
  if ( total_KE_CM <= 0. ) return 0.;
  const auto state = make_working_state( total_KE_CM, fragment_pdg,
    target_mass_for_charge(target_charge) );
  std::complex<double> S = s_matrix_element(fragment_pdg, two_j, l, two_s,
    state);

  // Guard against ±inf or NaN values that can occur in edge cases when the
  // Coulomb wavefunctions get huge, e.g., for low-energy alpha emission.
  // Numerical precision problems can lead to wrong answers, such as S == (inf,
  // 0) instead of the correct (1, 0).
  bool S_is_finite = std::isfinite( S.real() ) && std::isfinite( S.imag() );
  // If we have a ±inf or NaN in one of the components of S, then set the
  // transmission coefficient to zero
  if ( !S_is_finite ) return 0.;

  // To guard against numerical issues that can make the norm of the S-matrix
  // element creep above unity, explicitly enforce that it lies on the interval
  // [0, 1].
  // TODO: revisit this, perhaps add a warning message?
  double norm_S = std::norm(S);
  if ( norm_S < 0. || norm_S > 1.0000001 ) {
    MARLEY_LOG_DEBUG() << "Invalid S-matrix norm = " << norm_S << '\n';
  }
  norm_S = std::min(1., std::max(0., norm_S));

  // We can now compute the transmission coefficient in the usual way
  return 1.0 - norm_S;
}

std::complex<double>
marley::KoningDelarocheOpticalModel::s_matrix_element(int fragment_pdg,
  int two_j, int l, int two_s, WorkingState state) const
{
  // Compute request-local optical model parameters for the given fragment,
  // energy, and angular momenta.
  calculate_om_parameters(fragment_pdg, two_j, l, two_s, state);

  double step_size2_over_twelve = std::pow(step_size_, 2) / 12.0;

  std::complex<double> u1 = 0, u2 = 0;

  std::complex<double> a_n_minus_two;

  // a(r) really blows up at the origin for the optical model potential, but
  // we're saved by the boundary condition that u(0) = 0. We just need
  // something finite here, but we might as well make it zero.
  std::complex<double> a_n_minus_one = 0;
  std::complex<double> a_n = a(step_size_, l, state);

  std::complex<double> u_n_minus_two;
  // Boundary condition that the wavefunction vanishes at the origin (the
  // optical model potential blows up at r = 0)
  std::complex<double> u_n_minus_one = 0;

  // Asymptotic approximation for a regular potential (see J. Thijssen,
  // Computational Physics, p. 20 for details). We really just need something
  // finite and nonzero here, since our specific choice only determines the
  // overall normalization, which isn't important for determining the
  // transmission coefficients.
  std::complex<double> u_n = std::pow(step_size_, l + 1);

  // Optical model potential with and without the Coulomb potential included
  std::complex<double> U, U_minus_Vc;

  double r = step_size_;
  do {
    r += step_size_;
    a_n_minus_two = a_n_minus_one;
    a_n_minus_one = a_n;

    U_minus_Vc = omp_minus_Vc(r, state),
    U = U_minus_Vc + Vc(r, Rc, state.z, Z_);
    a_n = a(r, l, U, state);

    u_n_minus_two = u_n_minus_one;
    u_n_minus_one = u_n;

    u_n = ((2.0 - 10*step_size2_over_twelve*a_n_minus_one)*u_n_minus_one
      - (1.0 + step_size2_over_twelve*a_n_minus_two)*u_n_minus_two)
      / (1.0 + step_size2_over_twelve*a_n);
  }
  while (std::abs(U_minus_Vc) > MATCHING_RADIUS_THRESHOLD);

  double r_match_1 = r;
  u1 = u_n;

  /// @todo TODO: consider using a more sophisticated method for choosing the
  /// second matching radius
  // Advance at least as far as r_max. The actual maximum value used (which
  // will be an integer multiple of the step_size_) will be assigned to
  // r_match_2.
  double r_max = 1.2 * r_match_1;

  do {
    r += step_size_;
    a_n_minus_two = a_n_minus_one;
    a_n_minus_one = a_n;
    a_n = a(r, l, state);

    u_n_minus_two = u_n_minus_one;
    u_n_minus_one = u_n;

    u_n = ((2.0 - 10*step_size2_over_twelve*a_n_minus_one)*u_n_minus_one
      - (1.0 + step_size2_over_twelve*a_n_minus_two)*u_n_minus_two)
      / (1.0 + step_size2_over_twelve*a_n);
  }
  while (r < r_max);

  double r_match_2 = r;
  u2 = u_n;

  // Coulomb (Sommerfeld) parameter
  // Note that the relative (dimensionless) speed of the two particles
  // is just the speed of the fragment in the lab frame
  double beta_rel = marley_utils::real_sqrt( std::pow(state.fragment_KE_lab, 2)
    + 2.*state.fragment_KE_lab*state.fragment_mass )
    / (state.fragment_KE_lab + state.fragment_mass);

  // If beta_rel == 0, then eta blows up, so use a really small value
  /// @todo TODO: revisit this to see if you want to do something else
  if (beta_rel <= 0) beta_rel = 1e-8;

  double eta = Z_ * state.z * marley_utils::alpha / beta_rel;

  // Compute the Coulomb wavefunctions at the matching radii
  std::complex<double> Hplus1, Hminus1, Hplus2, Hminus2;

  // Fragment's CM frame wavenumber
  double k = marley_utils::real_sqrt( state.CM_frame_momentum_squared )
    / marley_utils::hbar_c;

  Hplus1 = coulomb_H_plus(l, eta, k*r_match_1);

  // H+ and H- are complex conjugates of each other
  Hminus1 = std::conj(Hplus1);

  Hplus2 = coulomb_H_plus(l, eta, k*r_match_2);
  Hminus2 = std::conj(Hplus2);

  // Compute the S matrix element using the radial wavefunction
  // evaluated at the two matching radii
  std::complex<double> S = (u1*Hminus2 - u2*Hminus1) / (u1*Hplus2 - u2*Hplus1);
  return S;
}

// Version of Schrodinger equation terms with the optical model potential
// U pre-computed
std::complex<double> marley::KoningDelarocheOpticalModel::a(double r,
  int l, std::complex<double> U, const WorkingState& state) const
{
  return (-l*(l+1) / std::pow(r, 2)) +
    (1. - (U / state.total_CM_frame_KE))
    * state.CM_frame_momentum_squared
    / marley_utils::hbar_c2;
}

// Non-derivative radial Schrödinger equation terms to use for computing
// transmission coefficients via the Numerov method
std::complex<double> marley::KoningDelarocheOpticalModel::a(double r, int l,
  const WorkingState& state) const
{
  return (-l*(l+1) / std::pow(r, 2)) +
    (1. - (omp(r, state) / state.total_CM_frame_KE))
    * state.CM_frame_momentum_squared
    / marley_utils::hbar_c2;
}

// Coulomb potential for a point particle with charge q*e interacting
// with a uniformly charged sphere with radius R and charge Q*e
double marley::KoningDelarocheOpticalModel::Vc(double r, double R, int Q,
  int q) const
{
  if (Q == 0 || q == 0) return 0.;
  else if (r < R) return Q * q * marley_utils::e2
    * (3. - std::pow(r / R, 2)) / (2. * R);
  else return Q * q * marley_utils::e2 / r;
}

// Woods-Saxon shape
double marley::KoningDelarocheOpticalModel::f(double r, double R, double a)
  const
{
  return std::pow(1 + std::exp((r - R) / a), -1);
}

// Compute the optical model potential at radius r
std::complex<double> marley::KoningDelarocheOpticalModel::omp(double r,
  const WorkingState& state) const
{
  return omp_minus_Vc(r, state) + Vc(r, Rc, state.z, Z_);
}

// Partial derivative with respect to r of the Woods-Saxon shape
double marley::KoningDelarocheOpticalModel::dfdr(double r, double R, double a) const
{
  // In the limit as r -> +-infinity, this goes to zero.
  // We pick an upper limit for the exponent to avoid evaluating
  // the function explicitly when r gets too large. Otherwise, C++
  // returns NaN because the function becomes indeterminate in double
  // precision (infinity/infinity or 0/0)
  double exponent = (r - R) / a;
  if (std::abs(exponent) > 100.) return 0;
  double temp = std::exp(exponent);
  return -temp / (a * std::pow(1 + temp, 2));
}

marley::KoningDelarocheOpticalModel::WorkingState
marley::KoningDelarocheOpticalModel::make_working_state(double KE_tot_CM,
  int fragment_pdg, double target_mass) const
{
  WorkingState state;
  state.target_mass = target_mass;
  calculate_kinematic_variables( KE_tot_CM, fragment_pdg, state );
  return state;
}

void marley::KoningDelarocheOpticalModel::calculate_kinematic_variables(
  double KE_tot_CM, int fragment_pdg, WorkingState& state) const
{
  // Store the total kinetic energy in the CM frame
  state.total_CM_frame_KE = KE_tot_CM;

  // Calculate the lab frame kinetic energy of the fragment from the
  // total CM frame kinetic energy
  const auto& mt = marley::MassTable::Instance();
  state.fragment_mass = mt.get_particle_mass( fragment_pdg );
  state.fragment_KE_lab = state.total_CM_frame_KE * (
    2.*(state.fragment_mass + state.target_mass) + state.total_CM_frame_KE )
    / (2. * state.target_mass);

  // Calculate the square of the CM frame 3-momentum of either particle
  state.CM_frame_momentum_squared = std::pow(state.target_mass, 2)
    * state.fragment_KE_lab * (2.*state.fragment_mass + state.fragment_KE_lab)
    / ( std::pow(state.fragment_mass + state.target_mass, 2)
    + 2.*state.target_mass*state.fragment_KE_lab );
}

double marley::KoningDelarocheOpticalModel::target_mass_for_charge(
  int target_charge) const
{
  // Compute the target mass based on its charge state without mutating the
  // shared optical-model object.
  const auto& mt = marley::MassTable::Instance();
  return mt.get_atomic_mass(Z_, A_)
    - target_charge*mt.get_particle_mass( marley_utils::ELECTRON );
}
