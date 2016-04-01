#include "meta_numerics.hh"
#include "Logger.hh"
#include "SphericalOpticalModel.hh"

std::complex<double> marley::SphericalOpticalModel::optical_model_potential(
  double r, double E, int fragment_pid, int two_j, int l, int two_s)
{
  calculate_om_parameters(E, fragment_pid, two_j, l, two_s);
  return omp(r);
}

// Finish an optical model potential calculation by taking the r
// dependence into account. Don't add in the Coulomb potential.
std::complex<double> marley::SphericalOpticalModel::omp_minus_Vc(double r) const
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
void marley::SphericalOpticalModel::calculate_om_parameters(double E,
  int fragment_pid, int two_j, int l, int two_s)
{
  // Fragment atomic, mass, and neutron numbers
  z = marley::MassTable::get_particle_Z(fragment_pid);
  int a = marley::MassTable::get_particle_A(fragment_pid);
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

marley::SphericalOpticalModel::SphericalOpticalModel(int z, int a) {
  Z = z;
  A = a;

  int N = A - Z; // Neutron number

  double A_to_the_one_third = std::pow(A, 1.0/3.0);

  // Initialize the spherical optical model parameters (see TALYS 1.6 manual)

  // Neutrons
  v1n = 59.30 - 21.0*(N - Z)/A - 0.024*A; // MeV
  v2n = 0.007228 - 1.48e-6*A; // MeV^(-1)
  v3n = 1.994e-5 - 2.0e-8*A; // MeV^(-2)
  v4n = 7e-9; // MeV^(-3)
  w1n = 12.195 + 0.0167*A; // MeV
  w2n = 73.55 + 0.0795*A; // MeV
  d1n = 16.0 - 16.0*(N - Z)/A; // MeV
  d2n = 0.0180 + 0.003802/(1 + std::exp((A - 156.)/8.0)); // MeV^(-1)
  d3n = 11.5; // MeV
  vso1n = 5.922 + 0.0030*A; // MeV
  vso2n = 0.0040; // MeV^(-1)
  wso1n = -3.1; // MeV
  wso2n = 160.; // MeV
  Efn = -11.2814 + 0.02646*A; // MeV
  Rvn = 1.3039*A_to_the_one_third - 0.4054; // fm
  avn = 0.6778 - 1.487e-4*A; // fm
  Rdn = 1.3424*A_to_the_one_third - 0.01585*std::pow(A_to_the_one_third, 2); // fm
  adn = 0.5446 - 1.656e-4*A; // fm
  Rso_n = 1.1854*A_to_the_one_third - 0.647; // fm
  aso_n = 0.59; // fm

  // Protons
  v1p = 59.30 + 21.0*(N - Z)/A - 0.024*A; // MeV
  v2p = 0.007067 + 4.23e-6*A; // MeV^(-1)
  v3p = 1.729e-5 + 1.136e-8*A; // MeV^(-2)
  v4p = v4n; // MeV^(-3)
  w1p = 14.667 + 0.009629*A; // MeV
  w2p = w2n; // MeV
  d1p = 16.0 + 16.0*(N - Z)/A; // MeV
  d2p = d2n; // MeV^(-1)
  d3p = d3n; // MeV
  vso1p = vso1n; // MeV
  vso2p = vso2n; // MeV^(-1)
  wso1p = wso1n; // MeV
  wso2p = wso2n; // MeV
  Efp = -8.4075 + 0.01378*A; // MeV
  Rvp = Rvn;
  avp = avn;
  Rdp = Rdn;
  adp = 0.5187 + 5.205e-4*A; // fm
  Rso_p = Rso_n;
  aso_p = aso_n;
  Rc = 1.198*A_to_the_one_third + 0.697/A_to_the_one_third
    + 12.994*std::pow(A_to_the_one_third, -4); // fm
  Vcbar_p = 1.73 * Z / Rc; // MeV

  // Compute and store reduced masses for all of the fragments of interest
  double m_nucleus = marley::MassTable::get_atomic_mass(Z, A);
  for (auto pid : marley::MassTable::get_fragment_pids()) {
    double m_fragment = marley::MassTable::get_particle_mass(pid);
    reduced_masses[pid] = (m_fragment * m_nucleus) / (m_fragment + m_nucleus);
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
double marley::SphericalOpticalModel::total_cross_section(double E,
  int fragment_pid, int two_s, size_t l_max, double h) //const
{
  double sum = 0.;
  for (size_t l = 0; l <= l_max; ++l) {
    int two_l = 2*l;
    for (int two_j = std::abs(two_l - two_s);
      two_j <= two_l + two_s; two_j += 2)
    {
      std::complex<double> S = s_matrix_element(E, fragment_pid, two_j, l,
        two_s, h);
      sum += (two_j + 1) * (1 - S.real());
    }
  }

  double mu = get_fragment_reduced_mass(fragment_pid);
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

double marley::SphericalOpticalModel::transmission_coefficient(double E,
  int fragment_pid, int two_j, int l, int two_s, double h) //const
{
  std::complex<double> S = s_matrix_element(E, fragment_pid, two_j, l,
    two_s, h);
  return 1.0 - std::norm(S);
}

std::complex<double> marley::SphericalOpticalModel::s_matrix_element(double E,
  int fragment_pid, int two_j, int l, int two_s, double h) //const
{

  // Update the optical model parameters stored in this object for the
  // given fragment, energy, and angular momenta
  calculate_om_parameters(E, fragment_pid, two_j, l, two_s);

  double h2_over_twelve = std::pow(h, 2) / 12.0;

  std::complex<double> u1 = 0, u2 = 0;

  std::complex<double> a_n_minus_two;
  // a(r) really blows up at the origin for the optical model potential, but we're saved
  // by the boundary condition that u(0) = 0. We just need something finite here, but we might
  // as well make it zero.
  std::complex<double> a_n_minus_one = 0;
  std::complex<double> a_n = a(h, E, fragment_pid, l);

  std::complex<double> u_n_minus_two;
  // Boundary condition that the wavefunction vanishes at the origin (the optical model
  // potential blows up at r = 0)
  std::complex<double> u_n_minus_one = 0;
  // Asymptotic approximation for a regular potential (see J. Thijssen, Computational
  // Physics, p. 20 for details). We really just need something finite and nonzero here, since
  // our specific choice only determines the overall normalization, which isn't important for
  // determining the transmission coefficients.
  std::complex<double> u_n = std::pow(h, l + 1);

  // Optical model potential with and without the Coulomb potential included
  std::complex<double> U, U_minus_Vc;

  double r = h;
  do {
    r += h;
    a_n_minus_two = a_n_minus_one;
    a_n_minus_one = a_n;

    U_minus_Vc = omp_minus_Vc(r),
    U = U_minus_Vc + Vc(r, Rc, z, Z);
    a_n = a(r, E, fragment_pid, l, U);

    u_n_minus_two = u_n_minus_one;
    u_n_minus_one = u_n;

    u_n = ((2.0 - 10*h2_over_twelve*a_n_minus_one)*u_n_minus_one
      - (1.0 + h2_over_twelve*a_n_minus_two)*u_n_minus_two)
      / (1.0 + h2_over_twelve*a_n);
  }
  while (std::abs(U_minus_Vc) > MATCHING_RADIUS_THRESHOLD);

  double r_match_1 = r;
  u1 = u_n;

  // TODO: consider using a more sophisticated method for choosing the second
  // matching radius
  // Advance at least as far as r_max. The actual maximum value used (which
  // will be an integer multiple of the step size h) will be assigned to
  // r_match_2.
  double r_max = 1.2 * r_match_1;

  do {
    r += h;
    a_n_minus_two = a_n_minus_one;
    a_n_minus_one = a_n;
    a_n = a(r, E, fragment_pid, l);

    u_n_minus_two = u_n_minus_one;
    u_n_minus_one = u_n;

    u_n = ((2.0 - 10*h2_over_twelve*a_n_minus_one)*u_n_minus_one
      - (1.0 + h2_over_twelve*a_n_minus_two)*u_n_minus_two)
      / (1.0 + h2_over_twelve*a_n);
  }
  while (r < r_max);

  double r_match_2 = r;
  u2 = u_n;

  LOG_DEBUG << "r1 = " << r_match_1 << ", r2 = " << r_match_2 << std::endl;

  // Coulomb parameter
  double mu = get_fragment_reduced_mass(fragment_pid);
  double k = marley_utils::real_sqrt(2.0 * mu * E) / marley_utils::hbar_c;
  // If k == 0, then eta blows up, so use a really small k instead of zero
  // (or a negative k due to roundoff error)
  if (k <= 0) k = 1e-8; //DEBUG!
  double eta = mu * Z * z * marley_utils::e2 / (marley_utils::hbar_c2 * k);

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
