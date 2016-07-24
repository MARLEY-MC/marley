#include "marley/Error.hh"
#include "marley/GammaStrengthFunctionModel.hh"
#include "marley/Level.hh"
#include "marley/Logger.hh"
#include "marley/Parity.hh"

using TrType = marley::GammaStrengthFunctionModel::TransitionType;

marley::GammaStrengthFunctionModel::GammaStrengthFunctionModel(int Z, int A)
  : Z_(Z), A_(A) {}

void marley::GammaStrengthFunctionModel::check_multipolarity(int l) {
  /// @todo Improve error message
  if (l < 1) throw marley::Error(std::string("Invalid multipolarity ")
    + std::to_string(l) + " given for gamma ray strength"
    + " function calculation");
}

TrType marley::GammaStrengthFunctionModel::determine_transition_type(
  int twoJi, marley::Parity Pi, marley::Level& level_f, int& l)
{
  int twoJf = level_f.twoJ();
  marley::Parity Pf = level_f.parity();

  return determine_transition_type(twoJi, Pi, twoJf, Pf, l);
}

TrType marley::GammaStrengthFunctionModel::determine_transition_type(
  marley::Level& level_i, marley::Level& level_f, int& l)
{
  int twoJi = level_i.twoJ();
  marley::Parity Pi = level_i.parity();

  int twoJf = level_f.twoJ();
  marley::Parity Pf = level_f.parity();

  return determine_transition_type(twoJi, Pi, twoJf, Pf, l);
}

// Returns whether a gamma transition from a nuclear state with spin twoJi / 2
// and parity Pi to a state with spin twoJf / 2 and parity Pf is electric or
// magnetic.  Also loads the integer l with the multipolarity of the
// transition.
TrType marley::GammaStrengthFunctionModel::determine_transition_type(
  int twoJi, marley::Parity Pi, int twoJf, marley::Parity Pf, int& l)
{
  if (twoJi == 0 && twoJf == 0) {
    MARLEY_LOG_WARNING() << "Unphysical Ji = 0 -> Jf = 0 EM transition"
      << " encountered in"
      << " GammaStrengthFunctionModel::determine_transition_type()";
    l = 1;
    return TrType::unphysical;
  }

  int two_delta_J = std::abs(twoJf - twoJi);

  // Odd values of two_delta_J are unphysical because they represent EM
  // transitions where the total angular momentum changes by half (photons are
  // spin 1)
  if (two_delta_J % 2) {
    MARLEY_LOG_WARNING() << "Unphysical EM transition between nuclear levels"
      << "with spins 2*Ji = " << std::to_string(twoJi) << " and 2*Jf = "
      << std::to_string(twoJf) << " encountered in"
      << " GammaStrengthFunctionModel::determine_transition_type()";
    l = two_delta_J / 2;
    return TrType::unphysical;
  }

  marley::Parity Pi_times_Pf = Pi * Pf;

  // Determine the multipolarity of this transition. Load l with the result.
  if (two_delta_J == 0) l = 1;
  // We already checked for unphysical odd values of two_delta_J above, so
  // using integer division here will always give the expected result.
  else l = two_delta_J / 2;

  // Determine whether this transition is electric or magnetic based on its
  // multipolarity and whether or not there is a change of parity

  // Pi * Pf = -1 gives electric transitions for odd l
  marley::Parity electric_parity;
  if (l % 2) electric_parity = -1; // l is odd
  // Pi * Pf = +1 gives electric transitions for even l
  else electric_parity = 1; // l is even

  if (Pi_times_Pf == electric_parity) return TrType::electric;
  else return TrType::magnetic;
}

double marley::GammaStrengthFunctionModel::transmission_coefficient(double Exi,
  int twoJi, marley::Parity Pi, marley::Level& level_f)
{
  double Exf = level_f.energy();
  int twoJf = level_f.twoJ();
  marley::Parity Pf = level_f.parity();

  // 0->0 EM transitions aren't allowed due to angular momentum conservation
  // (photons are spin 1), so skip the calculation if the initial and final
  // spins are zero.
  if (twoJi == 0 && twoJf == 0) return 0.;

  // If the final level has an energy greater than or equal to the initial
  // level, then this gamma-ray transition is trivially forbidden by energy
  // conservation.
  if (Exf >= Exi) return 0.;

  // Approximate the gamma energy by the energy difference between the two
  // levels
  // TODO: consider adding a nuclear recoil correction here
  double e_gamma = Exi - Exf;

  // Determine the type (electric or magnetic) and multipolarity of the gamma
  // transition between these two levels
  int mpol; // Multipolarity
  auto type = determine_transition_type(twoJi, Pi, twoJf, Pf, mpol);

  // The transmission coefficient below will be calculated differently by each
  // of the models that inherit from this class, since the version of
  // transmission_coefficient() used below is a virtual function.
  return transmission_coefficient(type, mpol, e_gamma);
}
