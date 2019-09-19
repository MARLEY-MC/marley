#include "marley/marley_kinematics.hh"
#include "marley/Reaction.hh"

using namespace marley_utils;
using ProcType = marley::Reaction::ProcessType;

namespace {

  std::map<ProcType, std::string> proc_type_to_string_map {
    { ProcType::NeutrinoCC, "\u03BD CC" },
    { ProcType::AntiNeutrinoCC, "anti-\u03BD CC" },
    { ProcType::NC, "NC" },
  };

  // Defines the neutrino species that can participate in each type
  // of scattering process
  std::map<ProcType, std::vector<int> > proc_type_to_nu_pdg = {

    { ProcType::NeutrinoCC,
      { ELECTRON_NEUTRINO, MUON_NEUTRINO, TAU_NEUTRINO }
    },

    { ProcType::AntiNeutrinoCC,
      { ELECTRON_ANTINEUTRINO, MUON_ANTINEUTRINO, TAU_ANTINEUTRINO }
    },

    { ProcType::NC,
      { ELECTRON_NEUTRINO, MUON_NEUTRINO, TAU_NEUTRINO,
        ELECTRON_ANTINEUTRINO, MUON_ANTINEUTRINO, TAU_ANTINEUTRINO }
    }

  };

} // Anonymous namespace


// Performs kinematics calculations for a two-two scattering reaction
// (a + b -> c + d)
void marley::Reaction::two_two_scatter(double KEa, double& s, double& Ec_cm,
  double& pc_cm, double& Ed_cm)
{
  // Get the lab-frame total energy of the projectile
  double Ea = KEa + ma_;

  // Compute Mandelstam s (the square of the total CM frame energy)
  s = ma_*ma_ + mb_*mb_ + 2.*mb_*Ea;
  double sqrt_s = std::sqrt(s);

  // Determine the CM frame energy and momentum of the ejectile
  Ec_cm = (s + mc_*mc_ - md_*md_) / (2 * sqrt_s);
  pc_cm = real_sqrt(std::pow(Ec_cm, 2) - mc_*mc_);

  // Determine the residue's CM frame energy. Roundoff errors may cause Ed_cm to
  // dip below md, which is unphysical. Prevent this from occurring by allowing
  // md to be the minimum value of Ed_cm. Also note that, in the CM frame, the
  // residue and ejectile have equal and opposite momenta.
  Ed_cm = std::max(sqrt_s - Ec_cm, md_);
}

marley::Event marley::Reaction::make_event_object(double KEa,
  double pc_cm, double cos_theta_c_cm, double phi_c_cm,
  double Ec_cm, double Ed_cm, double E_level)
{
  double sin_theta_c_cm = real_sqrt(1
    - std::pow(cos_theta_c_cm, 2));

  // Determine the Cartesian components of the ejectile's CM frame momentum
  double pc_cm_x = sin_theta_c_cm * std::cos(phi_c_cm) * pc_cm;
  double pc_cm_y = sin_theta_c_cm * std::sin(phi_c_cm) * pc_cm;
  double pc_cm_z = cos_theta_c_cm * pc_cm;

  // Get the lab-frame total energy of the projectile
  double Ea = KEa + ma_;

  // Determine the magnitude of the lab-frame 3-momentum of the projectile
  double pa = real_sqrt(KEa*(KEa + 2*ma_));

  // Create particle objects representing the projectile and target in the lab
  // frame
  // @todo Allow for projectile directions other than along the z-axis
  marley::Particle projectile(pdg_a_, Ea, 0, 0, pa, ma_);
  marley::Particle target(pdg_b_, mb_, 0, 0, 0, mb_);

  // Create particle objects representing the ejectile and residue in the CM
  // frame.
  marley::Particle ejectile(pdg_c_, Ec_cm, pc_cm_x, pc_cm_y, pc_cm_z, mc_);
  marley::Particle residue(pdg_d_, Ed_cm, -pc_cm_x, -pc_cm_y, -pc_cm_z, md_);

  // Boost the ejectile and residue into the lab frame.
  double beta_z = pa / (Ea + mb_);
  marley_kinematics::lorentz_boost(0, 0, -beta_z, ejectile);
  marley_kinematics::lorentz_boost(0, 0, -beta_z, residue);

  // Create the event object and load it with the appropriate information
  marley::Event event(projectile, target, ejectile, residue, E_level);
  return event;
}

int marley::Reaction::get_ejectile_pdg(int pdg_a, ProcType proc_type) {
  int pdg_c = 0;

  // First, check that the projectile PDG code is valid for the
  // given process type
  const auto& vec = proc_type_to_nu_pdg.at( proc_type );
  if ( std::find(vec.cbegin(), vec.cend(), pdg_a) != vec.end() ) {
    if ( proc_type == ProcType::NeutrinoCC ) pdg_c = pdg_a - 1;
    else if ( proc_type == ProcType::AntiNeutrinoCC ) pdg_c = pdg_a + 1;
    else if ( proc_type == ProcType::NC ) pdg_c = pdg_a;
    else throw marley::Error("Unrecognized ProcessType encountered in"
      " marley::Reaction::get_ejectile_pdg()");
  }
  else throw marley::Error("A projectile with PDG code "
    + std::to_string(pdg_a) + " cannot participate in reactions of type "
    + proc_type_to_string_map.at(proc_type));

  return pdg_c;
}

std::string marley::Reaction::proc_type_to_string(const ProcType& pt) {
  return proc_type_to_string_map.at( pt );
}

const std::vector<int>& marley::Reaction::get_projectiles(ProcType pt) {
  return proc_type_to_nu_pdg.at( pt );
}
