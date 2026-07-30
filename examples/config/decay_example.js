// Example MARLEY job configuration file for the "marley decay" command
// Steven Gardiner <gardiner@fnal.gov>
// Revised 19 July 2026 for MARLEY 2.0.0
//
// INTRODUCTION
//
// The "marley decay" command simulates stand-alone nuclear de-excitation
// cascades. Unlike "marley generate", it does NOT simulate a primary
// neutrino-nucleus scattering reaction. Instead, it starts from a
// user-specified nuclear excited state (defined by target nuclide,
// excitation energy, spin, and parity) and simulates the subsequent
// de-excitation cascade using MARLEY's nuclear de-excitation models.
//
// This is useful for studying the de-excitation behavior of specific
// nuclear levels in isolation, for generating de-excitation-only event
// samples, or for producing input to detector simulation frameworks.
//
// CONFIGURATION STRUCTURE
//
// A decay config file uses the same top-level format as the generate
// command and is parsed by the same JSONConfig class. The "reactions"
// key must be present but may be set to null (unlike for generate,
// where it must contain at least one reaction file). The "source" key
// may be omitted or set to null. The "generate" key is
// also ignored. Instead, a separate "decay" top-level object
// controls the run parameters for the decay command.
//
{ // An opening curly brace begins the configuration file content

  // RANDOM NUMBER SEED (optional)
  //
  // The "seed" key provides a 64-bit unsigned integer that will be used
  // to seed MARLEY's random number generator. If omitted, the system
  // time since the Unix epoch is used.
  seed: 123456,

  // REACTION INPUT FILES (required; may be null for decay)
  //
  // The "reactions" key must be present in the configuration file and
  // have a value that is either a JSON array of reaction file names
  // or null. Set it to null to skip reaction file loading, which is
  // the recommended setting for the decay command because no primary
  // neutrino-nucleus reactions are simulated.
  reactions: null,

  // NEUTRINO SOURCE SPECIFICATION (optional; may be null for decay)
  //
  // The "source" key may be omitted entirely or set to null. The decay
  // command does not rely on a neutrino source because it constructs
  // each excited nuclear state directly from the decays parameters below.
  source: null,

  // DECAY-SPECIFIC CONFIGURATION (required)
  //
  // The "decay" JSON object controls all aspects of the stand-alone
  // de-excitation simulation.
  //
  decay: {

    // EVENT COUNT (optional)
    //
    // The number of de-excitation events to generate before terminating.
    // Must be a positive integer. Default: 1000
    events: 10000,

    // PROJECTILE PDG CODE (optional)
    //
    // PDG code of the projectile that would cause the nuclear reaction
    // in a full scattering simulation. This is used only to determine
    // the ejectile PDG code via Reaction::get_ejectile_pdg(). The
    // resulting ejectile is stored in the HepMC3 event record for
    // consistency with the data format used by "marley generate".
    //
    // The default value is 12 (electron neutrino, nu_e).
    //
    // projectile: 12,

    // TARGET NUCLIDE (required)
    //
    // The target_Z and target_A keys define the nuclide that is
    // initially excited. For example, Z=18, A=40 corresponds to 40Ar.
    target_Z: 18,
    target_A: 40,

    // PROCESS TYPE (optional)
    //
    // An integer code representing the type of nuclear reaction that
    // populates the excited state. This determines the ejectile PDG
    // code and the change in proton number (Delta Z = projectile charge
    // - ejectile charge). The residue is defined as the daughter nucleus
    // after the ejectile has been emitted.
    //
    // Valid process types (only nuclear reaction types are allowed):
    //
    //   0  = NeutrinoCC_Discrete       1  = AntiNeutrinoCC_Discrete
    //   2  = NC_Discrete
    //   4  = NeutrinoCC_Continuum      5  = AntiNeutrinoCC_Continuum
    //   6  = NC_Continuum (default)
    //
    // The elastic scattering type (3 = NuElectronElastic) is not
    // allowed for the decay command.
    //
    // For CC reactions, the ejectile is a charged lepton (e-, mu-, or
    // tau-, depending on the projectile) and the residue has one more
    // proton than the target. For NC reactions, the ejectile is a
    // neutrino and the residue has the same Z as the target.
    //
    // proc_type: 6,

    // EXCITATION ENERGY (conditionally required)
    //
    // There are two ways to specify the excitation energy of the
    // initial nuclear state:
    //
    //   1. Fixed value: Use the "Ex" key with a single non-negative
    //      excitation energy in MeV.
    //
    //   2. Uniform sampling: Use both "Ex_min" and "Ex_max" keys.
    //      For each event, the excitation energy is sampled uniformly
    //      from the interval [Ex_min, Ex_max].
    //
    // Example: fixed excitation energy of 5.0 MeV
    Ex: 5.0,
    //
    // Alternative: uniform sampling from 4.0 to 6.0 MeV
    // Ex_min: 4.0, Ex_max: 6.0,

    // NUCLEAR SPIN (required)
    //
    // The "twoJ" key must be a JSON array of one or more positive
    // integers representing two times the total nuclear spin (so that
    // half-integer spins are represented). If multiple values are given,
    // the spin is sampled uniformly from the array for each event.
    //
    // Example: a single spin-1 state
    twoJ: [2],
    //
    // Example: uniform sampling of spin-0, spin-2, or spin-4
    // twoJ: [0, 4, 8],

    // INTRINSIC PARITY (optional)
    //
    // Parity of the initial nuclear state. Allowed values:
    //
    //   "+"      Positive parity (default)
    //   "-"      Negative parity
    //   "random" Randomly sample "+" or "-" with equal probability
    //            for each event
    //
    // parity: "+",

    // OUTPUT CONFIGURATION (optional)
    //
    // The output array follows the same format as the
    // generate.output array used by "marley generate".
    // Each entry is a JSON object with the following keys:
    //
    //   file    (string, optional)  Output file name.
    //                               Default: "decay_events.hepmc3"
    //
    //   format  (string, optional)  Output format. Valid values are
    //                               "ascii" and "root".
    //                               Default: "ascii"
    //
    //   mode    (string, optional)  File I/O mode. Valid values are
    //                               "overwrite" (default) and "resume".
    //                               Note: the decay command does not
    //                               implement resume/restore logic, so
    //                               the "resume" mode should not be used.
    //
    //   force   (bool,   optional)  Overwrite without prompting.
    //                               Default: false
    //
    output: [ { file: "decay_events.hepmc3", format: "ascii",
      mode: "overwrite" } ],
  },

  // Note on unused shared parameters:
  //
  // The following top-level keys accepted by JSONConfig are also
  // accepted by the decay command but have no effect on the
  // stand-alone de-excitation simulation: direction, target,
  // form_factors, generate, coulomb_mode,
  // do_deexcitations, sub_continuum_mode, energy_pdf_max, weights.
  //
  // Parameters that DO affect the de-excitation cascade:
  //
  //   opt_mod  (object, optional)
  //     Custom nuclear optical model parameters used in the
  //     Hauser-Feshbach decay width calculation. Alternative
  //     parameter sets will change the computed decay widths and
  //     branching ratios. See the "opt_mod" section of
  //     examples/config/annotated.js for details.
  //
  //   fragment_lmax  (int, optional, default 2)
  //     Maximum orbital angular momentum quantum number to consider
  //     when computing fragment (neutron, proton, alpha, etc.)
  //     decay widths in the continuum. Higher values increase
  //     computation time but may improve accuracy for high-energy
  //     transitions.
  //
  //   gamma_lmax  (int, optional, default 2, minimum 1)
  //     Maximum multipolarity to consider for gamma-ray decay
  //     widths. Higher values include higher-order electromagnetic
  //     transitions (E3, M3, etc.) at increased computational cost.

} // A closing curly brace should appear at the end of the file
