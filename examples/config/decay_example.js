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
// CONFIGURATION STRUCTURE
//
// A "marley decay" configuration file uses the same JSON-like format as
// the "marley generate" command (see examples/config/annotated.js for a
// full description of that format and its syntax rules). The "reactions"
// key must be present but may be set to null — unlike "marley generate",
// where it must contain at least one reaction input file name. The
// "source" key may be omitted entirely or set to null, and the "generate"
// key is ignored. A separate "decay" top-level object (described below)
// controls the run parameters specific to this command.
//
{ // An opening curly brace begins the configuration file content

  // RANDOM NUMBER SEED (optional)
  //
  // The "seed" key provides a nonnegative 63-bit integer (between 0
  // and 2^63 - 1, inclusive) that will be used to seed MARLEY's random
  // number generator.
  //
  // If this key is omitted, MARLEY will use the system time since the
  // Unix epoch as its random number seed.
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
  // each excited nuclear state directly from the decay parameters below.
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
    // The JSON parser expects this entry to be an integer literal, so
    // scientific notation is not currently allowed. If this key is
    // omitted, a value of 1000 will be assumed.
    events: 10000,

    // PROJECTILE PDG CODE (optional)
    //
    // PDG code of the notional projectile for this decay. Together with
    // the "proc_type" key (see below), this determines the PDG code of
    // the ejectile particle that is written into the output event record,
    // ensuring consistency with the HepMC3 event format produced by
    // "marley generate". For example, a projectile of 12 (νₑ) with a
    // CC process type yields an ejectile of 11 (e⁻), whereas a NC
    // process type yields an ejectile of 12 (νₑ). The projectile itself
    // is given zero momentum in the event record.
    //
    // The default value is 12 (electron neutrino, νₑ).
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
    // The "twoJ" key must be a JSON array of one or more nonnegative
    // integers representing two times the total nuclear spin (so that
    // half-integer spins can be represented as integers). If multiple
    // values are given, the spin is sampled uniformly from the array for
    // each event.
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
    // The "output" JSON array uses the same format as the
    // "generate.output" array described in examples/config/annotated.js,
    // with two differences: the default output file name is
    // "decay_events.hepmc3", and "resume" mode is not supported (an
    // error will be reported if it is requested).
    output: [ { file: "decay_events.hepmc3", format: "ascii",
      mode: "overwrite" } ],
  },

  // ADVANCED OPTIONS *********************************************************
  //
  // The following top-level keys from the "marley generate" configuration
  // are accepted by the decay command but have no effect on the stand-alone
  // de-excitation simulation: direction, target, form_factors, generate,
  // coulomb_mode, do_deexcitations, sub_continuum_mode, energy_pdf_max.
  // The status_update_interval key (available within the "generate" object
  // for "marley generate") is likewise not supported by "marley decay".
  //
  // The following top-level keys DO affect the de-excitation cascade:
  //
  // OPTICAL MODEL PARAMETERS (optional)
  //
  // The "opt_mod" key provides a custom configuration of nuclear optical
  // model parameters used in the Hauser-Feshbach decay width calculation.
  // Changing the parameter set alters the computed decay widths and
  // branching ratios. The format and available parameter sets are described
  // in the OPTICAL MODEL PARAMETERS section of examples/config/annotated.js
  // and in the full "opt_mod" documentation in
  // examples/config/reweight_example.js.
  //
  // opt_mod: #include:"optical_model_kduq_federal_cv.js"
  //
  // EVENT WEIGHT CALCULATORS (optional)
  //
  // The "weights" JSON array configures one or more weight calculators for the
  // de-excitation simulation, using the same types and format described in the
  // EVENT WEIGHT CALCULATORS section of examples/config/annotated.js. The
  // resulting weights are computed for each event and included in the output.
  //
  // ANGULAR MOMENTUM CUTOFFS (optional)
  //
  // The "fragment_lmax" key sets the maximum orbital angular momentum
  // quantum number to consider when computing fragment decay widths in the
  // continuum (default: 5). The "gamma_lmax" key sets the maximum
  // multipolarity for continuum gamma-ray decay widths (default: 5). Both
  // values are nonnegative integers; gamma_lmax must be >= 1.
  //
  // fragment_lmax: 5,
  // gamma_lmax: 5,

} // A closing curly brace should appear at the end of the file
