// Example MARLEY job configuration file for the "marley decay" command
// Steven Gardiner <gardiner@fnal.gov>
// Revised 20 July 2026 for MARLEY 2.0.0
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
// and "source" keys are not required and, if present, will be ignored.
// A separate "decay" top-level object (described below) controls the
// run parameters specific to this command.
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

    // INITIAL NUCLEAR STATE (required)
    //
    // The "nucleus" object specifies the nuclide and quantum numbers of
    // the excited nuclear state from which each de-excitation cascade
    // begins.
    nucleus: {

      // NUCLIDE SPECIFICATION (required)
      //
      // The nuclide may be specified using either (or both) of two
      // equivalent representations:
      //
      //   Option A: integer keys "Z" (proton number) and "A" (mass number)
      //   Option B: integer key "pdg" (nuclear PDG code = 10000*Z + 10*A
      //             + 1000000000; this is also the code stored in the
      //             output event record)
      //
      // At least one representation must be present. If both are given,
      // a marley::Error is thrown if they are inconsistent.
      //
      // Examples:
      //   Z: 18, A: 40          -> 40Ar  (nuclear PDG code 1000180400)
      //   pdg: 1000260560       -> 56Fe
      //   pdg: 1000060120       -> 12C
      //
      // Special cases: neutron (Z=0, A=1, pdg=2112) and proton
      //                (Z=1, A=1, pdg=2212) are also valid.
      //
      Z: 18,
      A: 40,
      // pdg: 1000180400,  // equivalent to Z=18, A=40

      // EXCITATION ENERGY (required; exactly one scheme)
      //
      // There are two ways to specify the excitation energy of the
      // initial nuclear state:
      //
      //   1. Fixed value: Use the "Ex" key with a single non-negative
      //      excitation energy in MeV.
      //
      //   2. Uniform sampling: Use both "Ex_min" and "Ex_max" keys.
      //      For each event the excitation energy is sampled uniformly
      //      from the interval [Ex_min, Ex_max] in MeV.
      //
      // When the excitation energy is below the unbound threshold for
      // the nuclide AND discrete level data are available in MARLEY's
      // structure database, the initial state is snapped to the nearest
      // tabulated discrete level. The level's spin and parity then
      // override the user-specified "twoJ" and "parity" values for that
      // event. A one-time warning is printed whenever a snap occurs.
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
      // values are given, the spin is sampled uniformly from the array
      // for each event.
      //
      // Note: when the initial state is snapped to a discrete level, the
      // level's spin overrides the sampled value (see EXCITATION ENERGY
      // above).
      //
      // Example: a single spin-1 state (2J = 2)
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
      // Note: when the initial state is snapped to a discrete level, the
      // level's parity overrides this value (see EXCITATION ENERGY above).
      //
      // parity: "+",

      // NET IONIC CHARGE (optional)
      //
      // Net ionic charge of the nucleus in units of the proton charge
      // (i.e., number of protons minus number of electrons). A value of
      // 0 corresponds to a neutral atom and is the default. A bare
      // (fully stripped) nucleus has net_charge = Z. A singly ionised
      // atom has net_charge = 1, and so on.
      //
      // The net charge is used to compute the nuclear mass from the
      // tabulated atomic mass: nuclear mass = atomic mass
      //   - net_charge * electron mass.
      //
      // Most users will want to leave this at the default.
      //
      // net_charge: 0,

    }, // end nucleus

    // OUTPUT CONFIGURATION (optional)
    //
    // The "output" JSON array uses the same format as the
    // "generate.output" array described in examples/config/annotated.js,
    // with two differences: the default output file name is
    // "decay_events.hepmc3", and "resume" mode is not supported (an
    // error will be reported if it is requested).
    output: [ { file: "decay_events.hepmc3", format: "ascii",
      mode: "overwrite" } ],

  }, // end decay

  // ADVANCED OPTIONS *********************************************************
  //
  // The following top-level keys from the "marley generate" configuration
  // are accepted by the decay command but have no effect on the stand-alone
  // de-excitation simulation: direction, target, form_factors, generate,
  // coulomb_mode, do_deexcitations, sub_continuum_mode, energy_pdf_max.
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
