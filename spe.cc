#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <eigen/Core>
#include <eigen/Eigenvalues>

#include "marley_utils.hh"
#include "meta_numerics.hh"
#include "TMarleyIntegrator.hh"
#include "TMarleyMassTable.hh"

class WoodsSaxon {
  public:
    inline WoodsSaxon(int z, int a) {
      if (z < 1) throw std::runtime_error(std::string("Cannot construct")
        + " a WoodsSaxon object for Z < 1");
      if (a < 1) throw std::runtime_error(std::string("Cannot construct")
        + " a WoodsSaxon object for A < 1");
      Z = z;
      A = a;
      double A_to_the_minus_one_third = std::pow(A, -1.0/3.0);
      hbar_omega = 45*A_to_the_minus_one_third
        - 25*std::pow(A_to_the_minus_one_third, 2);
      mN = 0.5 * (TMarleyMassTable::get_particle_mass(
        marley_utils::PROTON) + TMarleyMassTable::get_particle_mass(
        marley_utils::NEUTRON));
      b = marley_utils::hbar_c / std::sqrt(mN * hbar_omega);
      b2 = std::pow(b, 2);
      int N = A - Z;
      double factor = (static_cast<double>(N) - Z) / A;
      V0n = 51 - 33*factor; // MeV
      V0p = 51 + 33*factor; // MeV
      R = r0 / A_to_the_minus_one_third; // Nuclear radius (fm)
    }

    // The radial wavefunction for a 3D harmonic oscillator with principal
    // quantum number n (n >= 0), orbital angular momentum l (l >= 0), and
    // oscillator length b.
    double radial_3d_ho_wavefunction(int n, int l, double r) const;

    // Woods-Saxon shape
    inline double f(double r) const {
      return std::pow(1 + std::exp((r - R) / a), -1);
    }

    // Partial derivative with respect to r of the Woods-Saxon shape
    inline double dfdr(double r) const {
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

    // Coulomb potential for a point particle with charge z*e interacting with a
    // uniformly charged sphere with radius R and charge Z*e
    inline double Vc(double r) const {
      if (r < R) return Z * marley_utils::e2
        * (3 - std::pow(r / R, 2)) / (2 * R);
      else return Z * marley_utils::e2 / r;
    }

    inline double matrix_element(const TMarleyIntegrator& integrator,
      int two_j, int l, int pid, int nu_1, int nu_2) const
    {
      std::function<double(double)> f = std::bind(
        &WoodsSaxon::transformed_integrand, this, std::placeholders::_1,
        two_j, l, pid, nu_1, nu_2);
      return integrator.num_integrate(f, 0, 1);
    }

    inline double transformed_integrand(double t, int two_j, int l,
      int pid, int nu_1, int nu_2) const
    {
      if (t == 1 || t == 0) return 0.;
      else return matrix_element_integrand(t / (1 - t), two_j, l, pid, nu_1,
        nu_2) / std::pow(1 - t, 2);
    }

    inline double matrix_element_integrand(double r, int two_j, int l,
      int pid, int nu_1, int nu_2) const
    {
      double V0, vc;
      if (pid == marley_utils::NEUTRON) {
        V0 = V0n;
        vc = 0.;
      }
      else if (pid == marley_utils::PROTON) {
        V0 = V0p;
        vc = Vc(r);
      }
      else throw std::runtime_error(std::string("Woods-Saxon matrix elements")
        + " may only be computed for either protons or neutrons.");
      double v0ls = 0.44*V0;
      double kinetic_and_centrifugal_terms = (marley_utils::hbar_c2 / (2*mN))
        * ((4*nu_1 + 2*l + 3)/b2 - std::pow(r / b2, 2));
      double vws = -V0 * f(r);
      double vls = v0ls * (r02 / r) * dfdr(r);
      vls *= (two_j*(two_j + 2) - 4*l*(l + 1) - 3)/8.0;
      double g_nu_1_l = radial_3d_ho_wavefunction(nu_1, l, r);
      double g_nu_2_l = radial_3d_ho_wavefunction(nu_2, l, r);
      return std::pow(r, 2) * g_nu_1_l * g_nu_2_l
        * (kinetic_and_centrifugal_terms + vws + vc + vls);
    }

  private:
    int Z, A; // Atomic number and mass number
    double b; // Oscillator length (fm)
    double b2; // Oscillator length squared (fm^2)
    double hbar_omega; // Energy quantum (MeV)
    double mN; // Nucleon mass
    // Woods-Saxon parameters for the phenomenological nuclear potential
    constexpr static double r0 = 1.27; // fm
    constexpr static double r02 = std::pow(r0, 2); // fm^2
    constexpr static double a = 0.67; // Surface diffuseness (fm)
    double V0n, V0p; // Woods Saxon potential depths for neutrons and protons
    double R; // Nuclear radius
};

double WoodsSaxon::radial_3d_ho_wavefunction(int n, int l, double r) const {
  if (n < 0) throw std::runtime_error(std::string("Cannot compute")
    + " harmonic oscillator wavefunction for principal quantum number"
    + " n = " + std::to_string(n));
  if (l < 0) throw std::runtime_error(std::string("Cannot compute")
    + " harmonic oscillator wavefunction for orbital angular momentum"
    + " l = " + std::to_string(l));
  double num = meta_numerics::LogFactorial(n);
  double denom = meta_numerics::LogGamma(n + l + 1.5);
  double norm_factor = std::sqrt(2*std::exp(num - denom) / std::pow(b, 3));
  double r_over_b = r / b;
  double r_over_b_squared = std::pow(r_over_b, 2);
  return norm_factor * std::pow(r_over_b, l) * std::exp(-r_over_b_squared / 2)
    * meta_numerics::LaguerreL(n, l + 0.5, r_over_b_squared);
}

// Letters assigned to the values of the orbital angular momentum quantum
// number l (see https://en.wikipedia.org/wiki/Azimuthal_quantum_number for
// details)
const std::unordered_map<int, std::string> orbital_ang_momentum_letters = {
  { 0, "s" },
  { 1, "p" },
  { 2, "d" },
  { 3, "f" },
  { 4, "g" },
  { 5, "h" },
  { 6, "i" },
  { 7, "k" },
  { 8, "l" },
  { 9, "m" },
  { 10, "n" },
};

class SingleParticleEnergy {
  public:
    inline SingleParticleEnergy(double e, int nn, int twoj, int orbital_l,
      double sumA2)
    {
      E = e;
      n = nn;
      two_j = twoj; 
      l = orbital_l;
      sum_A2 = sumA2;
    }
    inline double get_energy() const { return E; }
    inline int get_n() const { return n; }
    inline int get_l() const { return l; }
    inline int get_two_j() const { return two_j; }
    inline double get_sum_A2() const { return sum_A2; }
  private:
    double E, sum_A2;
    int n, l, two_j;
};

int main() {

  std::unordered_map<int, std::vector<SingleParticleEnergy> > results;

  //std::cout << std::setprecision(16) << std::scientific;

  // Create a numerical integrator object to help compute the Hamiltonian
  // matrix elements.
  TMarleyIntegrator intg(50);

  constexpr int Z = 19;
  constexpr int A = 40;
  constexpr int l_max = 6;
  const std::vector<int> pids = { marley_utils::NEUTRON,
    marley_utils::PROTON };

  // Create the Woods-Saxon Hamiltonian object
  WoodsSaxon ws(Z, A);

  // max_nu*max_nu Hamiltonian matrix in the harmonic oscillator basis
  // Note that, since all of the matrix elements are real, we don't
  // need to worry about using a matrix type with std::complex elements
  constexpr size_t max_nu = 10;
  Eigen::MatrixXd hamiltonian_matrix(max_nu, max_nu);

  for (const auto pid : pids) {
    // Initialize the current list of results with an empty vector
    results[pid] = std::vector<SingleParticleEnergy>();
    for (int l = 0; l <= l_max; ++l) {
      for (int two_j = 2*l + 1; two_j >= std::abs(2*l - 1); two_j -= 2) {

        // Load the matrix with its elements for this iteration
        for (size_t i = 0; i < max_nu; ++i)
          for (size_t j = 0; j < max_nu; ++j)
            hamiltonian_matrix(i,j) = ws.matrix_element(intg, two_j, l, pid, i, j);

        // Solve for the eigensystem of the new Hamiltonian matrix
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eisolve(hamiltonian_matrix);
        if (eisolve.info() != Eigen::Success) {
          std::cout << "Fatal error: failed to diagonalize the Hamiltonian matrix."
            << std::endl;
          return 1;
        }

        //std::cout << "*** pid = " << pid << ", l = " << l << ", two_j = "
        //  << two_j << std::endl; std::cout << eisolve.eigenvalues() << std::endl;
        //std::cout << eisolve.eigenvectors() << std::endl;

        // If the diagonalization was successful, then print a report about the
        // single-particle energies.
        for (size_t c = 0; c < max_nu; ++c) {
          //if (l > static_cast<int>(c)) continue;
          double single_particle_energy = eisolve.eigenvalues()(c);
          double sum_squared_coeffs = 0.;
          for (size_t d = 0; d < max_nu; ++d)
            sum_squared_coeffs += std::pow(eisolve.eigenvectors().col(c)(d), 2);
          results.at(pid).push_back(SingleParticleEnergy(single_particle_energy,
            c, two_j, l, sum_squared_coeffs));
        }
      }
    }
  }

  // Print a report showing the computed single-particle energies
  for (auto& pair : results) {
    int pid = pair.first;
    std::string symbol = marley_utils::particle_symbols.at(pid);
    std::cout << "Table of single-particle energies for " << symbol
      << " states in " << A << marley_utils::element_symbols.at(Z)
      << std::endl;
    std::cout << "Orbital\t\tnl2j\t\t(n+1)l(j+1/2)\t\tEnergy (MeV)\t\tSum(A^2)"
      << std::endl;
    auto& vec = pair.second;
    // Sort the single-particle energies in ascending order
    std::sort(vec.begin(), vec.end(), [](const SingleParticleEnergy& spe1,
      const SingleParticleEnergy& spe2) -> bool { return spe1.get_energy()
      < spe2.get_energy(); });
    // Print the table
    for (const auto& spe : vec) {
      int n = spe.get_n();
      int two_j = spe.get_two_j();
      int l = spe.get_l();
      std::cout << n << orbital_ang_momentum_letters.at(l);
      if (two_j % 2) std::cout << two_j << "/2\t\t";
      else std::cout << two_j / 2 << "\t\t";
      std::cout << n << l << two_j << "\t\t"; 
      std::cout << n + 1 << l << static_cast<int>(
        std::round(two_j / 2.0 + 0.5)) << "\t\t";
      std::cout << spe.get_energy() << "\t\t" << spe.get_sum_A2() << std::endl;
    }
  }
  //std::cout << "The eigenvalues of H are:\n" << eisolve.eigenvalues()
  //  << std::endl;
  //std::cout << "Here's a matrix whose columns are eigenvectors of H \n"
  //  << "corresponding to these eigenvalues:\n"
  //  << eisolve.eigenvectors() << std::endl;
  //std::cout << eisolve.eigenvectors().col(0)(0) << std::endl;
  return 0;
}

  //std::function<double(double)> f = [&ws](double t) -> double {
  //  if (t == 1) return 0.;
  //  double r = t / (1 - t);
  //  double g1 = ws.radial_3d_ho_wavefunction(1, 0, r);
  //  double g2 = ws.radial_3d_ho_wavefunction(1, 0, r);
  //  return std::pow(r, 2) * g1 * g2 / std::pow(1 - t, 2);
  //};
  //std::cout << integrator.num_integrate(f, 0, 1) << std::endl;
