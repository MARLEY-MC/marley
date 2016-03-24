#include "DecayChannel.hh"
#include "NuclearPhysics.hh"

void marley::ContinuumFragmentDecayChannel::get_post_decay_parameters(double& Ex,
  int& two_J, marley::Parity& Pi)
{
  // Store the initial excitation energy for later use
  double Exi = Ex;

  // Sample a final excitation energy uniformly on [Emin, Emax)
  Ex = gen.uniform_random_double(Emin, Emax, include_upper_edge);

  // Sample a final spin and parity
  // Compute outgoing fragment kinetic energy in the rest frame of the
  // initial nucleus
  double Ea = (Mconst - Ex*(2*Mfgs + Ex)) / (2 * (Migs + Exi));

  marley::NuclearPhysics::sample_fragment_spin_parity(two_J, Pi, fragment,
    om, gen, Ex, Ea);
}

void marley::ContinuumGammaDecayChannel::get_post_decay_parameters(double& Ex,
  int& twoJ, marley::Parity& Pi)
{
  // Store the initial excitation energy for later use
  double Exi = Ex;

  // Sample a final excitation energy uniformly on [Emin, Emax)
  Ex = gen.uniform_random_double(Emin, Emax, include_upper_edge);

  // Sample a final spin and parity
  marley::NuclearPhysics::sample_gamma_spin_parity(Zi, Ai, twoJ, Pi, Exi,
    Ex, gen);
}
