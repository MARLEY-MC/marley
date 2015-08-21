#include "TMarleyNuclearPhysics.hh"

using TrType = TMarleyNuclearPhysics::TransitionType;

TrType TMarleyNuclearPhysics::determine_gamma_transition_type(
  TMarleyLevel* level_i, TMarleyLevel* level_f, int& l)
{
  int Ji = level_i->get_two_J() / 2;
  TMarleyParity Pi = level_i->get_parity();

  int Jf = level_f->get_two_J() / 2;
  TMarleyParity Pf = level_f->get_parity();

  return determine_gamma_transition_type(Ji, Pi, Jf, Pf, l);
}


// Returns whether a gamma transition from a nuclear state with spin Ji and
// parity Pi to a state with spin Jf and parity Pf is electric or magnetic.
// Also loads the integer l with the multipolarity of the transition.
TrType TMarleyNuclearPhysics::determine_gamma_transition_type(int Ji,
  int Pi, int Jf, int Pf, int& l)
{
  // TODO: reconsider how to handle this situation
  if (Ji == 0 && Jf == 0) throw std::runtime_error(
    std::string("0 -> 0 EM transitions are not allowed."));

  int delta_J = std::abs(Jf - Ji);
  int Pi_times_Pf = Pi * Pf;
  // TODO: add error check to verify that each parity is either +1 or -1

  // Determine the multipolarity of this transition. Load l with the result.
  if (delta_J == 0) l = 1;
  else l = delta_J;

  // Determine whether this transition is electric or magnetic based on its
  // multipolarity and whether or not there is a change of parity

  // Pi * Pf = -1 gives electric transitions for odd l
  int electric_parity;
  if (l % 2) electric_parity = -1; // l is odd
  // Pi * Pf = +1 gives electric transitions for even l
  else electric_parity = 1; // l is even

  if (Pi_times_Pf == electric_parity) return TransitionType::electric;
  else return TransitionType::magnetic;
}

double TMarleyNuclearPhysics::gamma_strength_function(int Z, int A,
  TransitionType type, int l, double e_gamma)
{
  // TODO: improve this error message
  if (l < 1) throw std::runtime_error(std::string("Invalid multipolarity")
    + std::to_string(l) + " given for gamma ray strength"
    + " function calculation");

  // The strength, energy, and width of the giant resonance for a transition of
  // type x (E or M) and multipolarity l
  double sigma_xl, e_xl, gamma_xl;

  if (type == TransitionType::electric) {
    if (l == 1) {
      e_xl = 31.2*std::pow(A, -1.0/3.0) + 20.6*std::pow(A, -1.0/6.0);
      gamma_xl = 0.026 * std::pow(e_xl, 1.91);
      sigma_xl = 1.2 * 120 * (A - Z) * Z/(A * marley_utils::pi * gamma_xl)
        * marley_utils::mb;
    }
    if (l > 1) {
      // Values for E2 transitions
      e_xl = 63*std::pow(A, -1.0/3.0);
      gamma_xl = 6.11 - 0.012*A;
      sigma_xl = 0.00014 * std::pow(Z, 2) * e_xl
        / (std::pow(A, 1.0/3.0) * gamma_xl) * marley_utils::mb;
      // If this is an E2 transition, we're done. Otherwise,
      // compute the giant resonance strength iteratively
      for (int i = 2; i < l; ++i) {
        sigma_xl *= 8e-4;
      }
    }
  }
  else if (type == TransitionType::magnetic) {
    // Values for M1 transitions
    // The commented-out version is for RIPL-1
    //double factor_m1 = 1.58e-9*std::pow(A, 0.47);
    // RIPL-2 factor
    const double e_gamma_ref = 7.0; // MeV
    double factor_m1 = gamma_strength_function(Z, A,
      TransitionType::electric, 1, e_gamma_ref)
      / (0.0588 * std::pow(A, 0.878));
    gamma_xl = 4.0;
    e_xl = 41*std::pow(A, -1.0/3.0);
    sigma_xl = (std::pow(std::pow(e_gamma_ref, 2) - std::pow(e_xl, 2), 2)
      + std::pow(e_gamma_ref, 2) * std::pow(gamma_xl, 2))
      * (3 * std::pow(marley_utils::pi, 2) * factor_m1)
      / (e_gamma_ref * std::pow(gamma_xl, 2));
    // If this is an M1 transition, we're done. Otherwise,
    // compute the giant resonance strength iteratively
    for (int i = 1; i < l; ++i) {
      sigma_xl *= 8e-4;
    }
  }
  // TODO: improve this error message
  else throw std::runtime_error(std::string("Invalid transition type")
    + " given for gamma ray strength function calculation");

  // Now that we have the appropriate giant resonance parameters,
  // calculate the strength function using the Brink-Axel expression.
  // Note that the strength function has units of MeV^(-3)
  double f_xl = (sigma_xl * std::pow(e_gamma, 3 - 2*l)
    * std::pow(gamma_xl, 2)) / ((2*l + 1) * std::pow(marley_utils::pi, 2)
    * (std::pow(std::pow(e_gamma, 2) - std::pow(e_xl, 2), 2)
    + std::pow(e_gamma, 2) * std::pow(gamma_xl, 2)));

  return f_xl;
}

// Computes the partial decay width for a gamma transition under the Weisskopf
// single-particle approximation.
double TMarleyNuclearPhysics::weisskopf_partial_decay_width(int A,
  TransitionType type, int l, double e_gamma)
{
  // Compute double factorial of 2l + 1
  int dfact = 1;
  for (int n = 2*l + 1; n > 0; n -= 2) dfact *= n;

  // Multipolarity factor
  double lambda = (l + 1) / (l * std::pow(dfact, 2))
    * std::pow(3.0 / (l + 3), 2);

  // Estimated nuclear radius (fm)
  double R = 1.2 * std::pow(A, 1.0/3.0);

  // Electric transition partial decay width
  double el_width = 2 * marley_utils::alpha * lambda
    * std::pow(R * e_gamma / marley_utils::hbar_c, 2*l) * e_gamma;

  if (type == TransitionType::electric) {
    return el_width;
  }

  else if (type == TransitionType::magnetic) {
    double mp = TMarleyMassTable::get_particle_mass(PROTON);
    return 10 * el_width * std::pow(marley_utils::hbar_c / (mp * R), 2);
  }

  // TODO: improve this error message
  else throw std::runtime_error(std::string("Invalid transition type")
    + " given for Weisskopf gamma-ray partial width calculation");
}
