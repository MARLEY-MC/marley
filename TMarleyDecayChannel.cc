#include "TMarleyDecayChannel.hh"
#include "TMarleyNuclearPhysics.hh"

void TMarleyContinuumFragmentDecayChannel::get_post_decay_parameters(double& Ex,
  int& two_J, TMarleyParity& Pi)
{
  // Store the initial excitation energy for later use
  double Exi = Ex;

  // Sample a final excitation energy uniformly on [Emin, Emax)
  Ex = gen.uniform_random_double(Emin, Emax, include_upper_edge);

  // Sample a final spin and parity
  // Compute outgoing fragment kinetic energy in the rest frame of the
  // initial nucleus
  double Ea = (Mconst - Ex*(2*Mfgs + Ex)) / (2 * (Migs + Exi));

  TMarleyNuclearPhysics::sample_fragment_spin_parity(two_J, Pi, fragment,
    om, gen, Ex, Ea);
}

void TMarleyContinuumGammaDecayChannel::get_post_decay_parameters(double& Ex,
  int& twoJ, TMarleyParity& Pi)
{
  // Store the initial excitation energy for later use
  double Exi = Ex;

  // Sample a final excitation energy uniformly on [Emin, Emax)
  Ex = gen.uniform_random_double(Emin, Emax, include_upper_edge);

  // Sample a final spin and parity
  TMarleyNuclearPhysics::sample_gamma_spin_parity(Zi, Ai, twoJ, Pi, Exi,
    Ex, gen);
}
