#include "marley/marley_utils.hh"
#include "marley/meta_numerics.hh"
#include "marley/HauserFeshbachDecay.hh"
#include "marley/Logger.hh"
#include "marley/KoningDelarocheOpticalModel.hh"

std::complex<double>
marley::KoningDelarocheOpticalModel::optical_model_potential(double r,
  double E, int fragment_pdg, int two_j, int l, int two_s)
{
  calculate_om_parameters(E, fragment_pdg, two_j, l, two_s);
  return omp(r);
}

// Finish an optical model potential calculation by taking the r
// dependence into account. Don't add in the Coulomb potential.
std::complex<double>
marley::KoningDelarocheOpticalModel::omp_minus_Vc(double r) const
{
  double f_v = f(r, Rv, av);
  double dfdr_d = dfdr(r, Rd, ad);

  double temp_Vv = Vv * f_v;
  double temp_Wv = Wv * f_v;
  double temp_Wd = -4 * Wd * ad * dfdr_d;

  double temp_Vso = 0;
  double temp_Wso = 0;

  if (spin_orbit_eigenvalue != 0) {

    double factor_so = lambda_piplus2 * dfdr(r, Rso, aso)
      * spin_orbit_eigenvalue / r;

    temp_Vso = Vso * factor_so;
    temp_Wso = Wso * factor_so;
  }

  return std::complex<double>(-temp_Vv + temp_Vso,
    -temp_Wv - temp_Wd + temp_Wso);
}

// Compute all of the pieces of the optical model that depend on the
// fragment's energy E but not on its distance from the origin r.
// Store them in the appropriate class members.
void marley::KoningDelarocheOpticalModel::calculate_om_parameters(double E,
  int fragment_pdg, int two_j, int l, int two_s)
{
  // Fragment atomic, mass, and neutron numbers
  z = marley_utils::get_particle_Z(fragment_pdg);
  int a = marley_utils::get_particle_A(fragment_pdg);
  int n = a - z;

  // Eigenvalue of the spin-orbit operator
  // 2*(l.s) = j*(j + 1)  - l*(l + 1) -  s*(s + 1)
  // = 0.25*((2j - 2s)*(2j + 2s + 2)) - l*(l+1)
  // (to keep the units right we take hbar = 1).
  bool spin_zero = two_s == 0;
  if (spin_zero) spin_orbit_eigenvalue = 0;
  else spin_orbit_eigenvalue = 0.25*((two_j - two_s)
    * (two_j + two_s + 2)) - l*(l + 1);


  // Geometrical parameters
  Rv = 0;
  av = 0;
  Rd = 0;
  ad = 0;
  Rso = 0;
  aso = 0;

  // Terms in the spherical optical model potential
  Vv = 0;
  Wv = 0;
  Wd = 0;
  Vso = 0;
  Wso = 0;

  // Energy to use when computing folded potentials
  double E_eff = E / a;

  if (n > 0) {
    double Ediff_n = E_eff - Efn;
    double Ediff_n2 = std::pow(Ediff_n, 2);
    double Ediff_n3 = std::pow(Ediff_n, 3);

    Vv += n * v1n * (1 - v2n*Ediff_n + v3n*Ediff_n2 - v4n*Ediff_n3);
    Wv += n * w1n * Ediff_n2 / (Ediff_n2 + std::pow(w2n, 2));
    Wd += n * d1n * Ediff_n2 * std::exp(-d2n * Ediff_n)
      / (Ediff_n2 + std::pow(d3n, 2));

    Rv += n * Rvn;
    av += n * avn;
    Rd += n * Rdn;
    ad += n * adn;
    Rso += n * Rso_n;
    aso += n * aso_n;

    if (!spin_zero) {
      double Ediff_so_n = E - Efn;
      double Ediff_so_n2 = std::pow(Ediff_so_n, 2);
      Vso += vso1n * std::exp(-vso2n * Ediff_so_n);
      Wso += wso1n * Ediff_so_n2 / (Ediff_so_n2 + std::pow(wso2n, 2));
    }
  }

  if (z > 0) {
    double Ediff_p = E_eff - Efp;
    double Ediff_p2 = std::pow(Ediff_p, 2);
    double Ediff_p3 = std::pow(Ediff_p, 3);

    Vv += z * v1p * (1 - v2p*Ediff_p + v3p*Ediff_p2 - v4p*Ediff_p3
      + Vcbar_p*(v2p - 2*v3p*Ediff_p + 3*v4p*Ediff_p2));
    Wv += z * w1p * Ediff_p2 / (Ediff_p2 + std::pow(w2p, 2));
    Wd += z * d1p * Ediff_p2 * std::exp(-d2p * Ediff_p)
      / (Ediff_p2 + std::pow(d3p, 2));

    Rv += z * Rvp;
    av += z * avp;
    Rd += z * Rdp;
    ad += z * adp;
    Rso += z * Rso_p;
    aso += z * aso_p;

    if (!spin_zero) {
      double Ediff_so_p = E - Efp;
      double Ediff_so_p2 = std::pow(Ediff_so_p, 2);
      Vso += vso1p * std::exp(-vso2p * Ediff_so_p);
      Wso += wso1p * Ediff_so_p2 / (Ediff_so_p2 + std::pow(wso2p, 2));
    }
  }

  if (a > 1) {
    Rv /= a;
    av /= a;
    Rd /= a;
    ad /= a;
    Rso /= a;
    aso /= a;

    // Apply folding factor for composite particle spin-orbit potentials
    if (!spin_zero) {
      bool z_odd = z % 2;
      bool n_odd = n % 2;
      // This factor stays zero for even-even nuclides (which should all be
      // spin-zero anyway)
      double factor = 0.;
      if (z_odd && n_odd) factor = 2.0; // odd-odd
      else if (z_odd != n_odd) factor = 1.0; // even-odd
      factor /= 2*a;

      Vso *= factor;
      Wso *= factor;
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

  // Compute and store reduced masses for all of the fragments of interest
  const marley::MassTable& mt = marley::MassTable::Instance();
  double m_nucleus = mt.get_atomic_mass(Z_, A_);
  for (const auto& f : marley::HauserFeshbachDecay::get_fragments()) {
    int pdg = f.get_pid();
    double m_fragment = mt.get_particle_mass(pdg);
    reduced_masses_[pdg] = (m_fragment * m_nucleus) / (m_fragment + m_nucleus);
  }
}

// Computes the total cross section (in barns) for a given fragment
// scattering on this nucleus at a given kinetic energy E.
//
// The total cross section given here by the optical model may be directly
// compared to experiment. Expressions exist for the optical model elastic and
// reaction (absorption) cross sections, but these are hard to directly compare
// to the data because the absorption cross section includes the compound
// elastic channel.
double marley::KoningDelarocheOpticalModel::total_cross_section(double E,
  int fragment_pdg, int two_s, size_t l_max)
{
  double sum = 0.;
  for (size_t l = 0; l <= l_max; ++l) {
    int two_l = 2*l;
    for (int two_j = std::abs(two_l - two_s);
      two_j <= two_l + two_s; two_j += 2)
    {
      std::complex<double> S = s_matrix_element(E, fragment_pdg, two_j, l,
        two_s);
      sum += (two_j + 1) * (1 - S.real());
    }
  }

  double mu = get_fragment_reduced_mass(fragment_pdg);
  double k = marley_utils::real_sqrt(2.0 * mu * E) / marley_utils::hbar_c;
  // If k == 0, then eta blows up, so use a really small k instead of zero
  // (or a negative k due to roundoff error)
  if (k <= 0) k = 1e-8; //DEBUG!

  // This expression is based on equation 3.2.34 from I. Thompson and F. Nunes,
  // Nuclear Reactions for Astrophysics, p. 85. Note that we have reduced the
  // coupled-channels form shown therein to a single channel, and we have also
  // assumed that the nucleus has zero spin (spherical optical model).
  //
  // Compute the cross section in units of fm^2
  double xs = marley_utils::two_pi * sum / ((two_s + 1) * std::pow(k, 2));
  // Return the cross section in barns (1 b = 100 fm^2)
  return xs / 100.;
}

double marley::KoningDelarocheOpticalModel::transmission_coefficient(double E,
  int fragment_pdg, int two_j, int l, int two_s)
{
  std::complex<double> S = s_matrix_element(E, fragment_pdg, two_j, l,
    two_s);
  return 1.0 - std::norm(S);
}

std::complex<double>
marley::KoningDelarocheOpticalModel::s_matrix_element(double E,
  int fragment_pdg, int two_j, int l, int two_s)
{

  // Update the optical model parameters stored in this object for the
  // given fragment, energy, and angular momenta
  calculate_om_parameters(E, fragment_pdg, two_j, l, two_s);

  double step_size2_over_twelve = std::pow(step_size_, 2) / 12.0;

  std::complex<double> u1 = 0, u2 = 0;

  std::complex<double> a_n_minus_two;

  // a(r) really blows up at the origin for the optical model potential, but
  // we're saved by the boundary condition that u(0) = 0. We just need
  // something finite here, but we might as well make it zero.
  std::complex<double> a_n_minus_one = 0;
  std::complex<double> a_n = a(step_size_, E, fragment_pdg, l);

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

    U_minus_Vc = omp_minus_Vc(r),
    U = U_minus_Vc + Vc(r, Rc, z, Z_);
    a_n = a(r, E, fragment_pdg, l, U);

    u_n_minus_two = u_n_minus_one;
    u_n_minus_one = u_n;

    u_n = ((2.0 - 10*step_size2_over_twelve*a_n_minus_one)*u_n_minus_one
      - (1.0 + step_size2_over_twelve*a_n_minus_two)*u_n_minus_two)
      / (1.0 + step_size2_over_twelve*a_n);
  }
  while (std::abs(U_minus_Vc) > MATCHING_RADIUS_THRESHOLD);

  double r_match_1 = r;
  u1 = u_n;

  // TODO: consider using a more sophisticated method for choosing the second
  // matching radius
  // Advance at least as far as r_max. The actual maximum value used (which
  // will be an integer multiple of the step_size_) will be assigned to
  // r_match_2.
  double r_max = 1.2 * r_match_1;

  do {
    r += step_size_;
    a_n_minus_two = a_n_minus_one;
    a_n_minus_one = a_n;
    a_n = a(r, E, fragment_pdg, l);

    u_n_minus_two = u_n_minus_one;
    u_n_minus_one = u_n;

    u_n = ((2.0 - 10*step_size2_over_twelve*a_n_minus_one)*u_n_minus_one
      - (1.0 + step_size2_over_twelve*a_n_minus_two)*u_n_minus_two)
      / (1.0 + step_size2_over_twelve*a_n);
  }
  while (r < r_max);

  double r_match_2 = r;
  u2 = u_n;

  // Coulomb parameter
  double mu = get_fragment_reduced_mass(fragment_pdg);
  double k = marley_utils::real_sqrt(2.0 * mu * E) / marley_utils::hbar_c;

  // If k == 0, then eta blows up, so use a really small k instead of zero
  // (or a negative k due to roundoff error)
  if (k <= 0) k = 1e-8; //DEBUG!

  double eta = mu * Z_ * z * marley_utils::e2 / (marley_utils::hbar_c2 * k);

  // Compute the Coulomb wavefunctions at the matching radii
  double F, G;
  std::complex<double> Hplus1, Hminus1, Hplus2, Hminus2;
  F = meta_numerics::CoulombF(l, eta, k*r_match_1);
  G = meta_numerics::CoulombG(l, eta, k*r_match_1);
  Hplus1 = std::complex<double>(G, F);

  // H+ and H- are complex conjugates of each other
  Hminus1 = std::conj(Hplus1);

  F = meta_numerics::CoulombF(l, eta, k*r_match_2);
  G = meta_numerics::CoulombG(l, eta, k*r_match_2);
  Hplus2 = std::complex<double>(G, F);
  Hminus2 = std::conj(Hplus2);

  // Compute the S matrix element using the radial wavefunction
  // evaluated at the two matching radii
  std::complex<double> S = (u1*Hminus2 - u2*Hminus1) / (u1*Hplus2 - u2*Hplus1);
  return S;
}

// Version of Schrodinger equation terms with the optical model potential
// U pre-computed
std::complex<double> marley::KoningDelarocheOpticalModel::a(double r,
  double E, int fragment_pdg, int l, std::complex<double> U) const
{
  return (-l*(l+1) / std::pow(r, 2)) +
    2 * reduced_masses_.at(fragment_pdg) * (E - U) / marley_utils::hbar_c2;
}

// Non-derivative radial Schrödinger equation terms to use for computing
// transmission coefficients via the Numerov method
std::complex<double> marley::KoningDelarocheOpticalModel::a(double r,
  double E, int fragment_pdg, int l)
{
  return (-l*(l+1) / std::pow(r, 2)) +
    2 * reduced_masses_.at(fragment_pdg) * (E - omp(r))
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
std::complex<double> marley::KoningDelarocheOpticalModel::omp(double r) const
{
  return omp_minus_Vc(r) + Vc(r, Rc, z, Z_);
}

// Partial derivative with respect to r of the Woods-Saxon shape
double marley::KoningDelarocheOpticalModel::dfdr(double r, double R, double a) const
{
  // In the limit as r -> +-infinity, this goes to zero.
  // We pick an upper limit for the exponent to avoid evaluating
  // the function explicitly when r gets too large (otherwise, C++
  // returns NaN because the function becomes indeterminate in double
  // precision (infinity/infinity or 0/0)
  double exponent = (r - R) / a;
  if (std::abs(exponent) > 100.) return 0;
  double temp = std::exp(exponent);
  return -temp / (a * std::pow(1 + temp, 2));
}
