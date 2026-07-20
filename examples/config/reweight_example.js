// Example MARLEY job configuration file for the "marley reweight" command
// Steven Gardiner <gardiner@fnal.gov>
// Revised 19 July 2026 for MARLEY 2.0.0
//
// INTRODUCTION
//
// The "marley reweight" command assigns weights to a previously generated
// sample of Monte Carlo events. It reads a configuration file that specifies
// one or more weight calculators, applies them to each event in one or more
// input HepMC3-format files, and writes the weights together with the original
// event data to an output file.
//
// The configuration format for the "marley reweight" command uses the same
// JSON-like syntax as the other MARLEY commands. The "weights" top-level key
// specifies the weight calculator(s) to run, while the optional "reweight"
// section may be used to configure output settings. All other parameters needed
// by the Generator are automatically restored from the "MARLEY.JSONconfig"
// attribute saved in the first input event file. For multiple input files,
// compatibility of their saved run configurations is checked, and the reweight
// job will abort if discrepancies are found.
//
{ // An opening curly brace begins the configuration file content

  // EVENT WEIGHT CALCULATORS (required)
  //
  // The "weights" JSON array configures one or more weight calculators that
  // assign additional event weights during reweighting. Each entry is a JSON
  // object with at least a "type" key and a "name" key. The three available
  // types are "trivial", "optical_model", and "strength_variation". Basic
  // documentation of the weight calculator types and their parameters is
  // provided in the annotated example configuration file for "marley generate"
  // (examples/config/annotated.js). The notes below provide some additional
  // details.
  weights: [

    // Produces a unit weight for every event
    { type: "trivial", name: "MyTrivialWeight" },

    // OPTICAL MODEL WEIGHT CALCULATOR
    //
    // An "optical_model" weight calculator recomputes the Hauser-Feshbach
    // decay widths for each event using a user-specified alternative set of
    // optical model parameters. Weights are evaluated based on ratios of
    // the alternative to the original decay widths. For events involving
    // multiple Hauser-Feshbach decay vertices, the per-vertex weights
    // are multiplied together to obtain an overall event weight.
    //
    // The required "opt_mod" key provides the alternative optical model
    // parameters. Its value is a JSON object containing one or more
    // named sub-objects. At minimum, a "Default" key must be present,
    // whose value supplies the 46 double-valued parameters of the
    // Koning-Delaroche optical model potential
    // (https://doi.org/10.1016/S0375-9474(02)01321-0). All energies below
    // are in MeV and all lengths are in fm. An optional "step_size" key
    // (in fm) controls the step size used for numerical solution of the
    // Schrödinger equation via the Numerov method.
    //
    //   Real volume potential (Vv):
    //     v10      (MeV)    -- v1 constant term
    //     v1A      (MeV)    -- v1 mass-number coefficient
    //     v1alpha  (MeV)    -- v1 asymmetry coefficient
    //     vn20     (MeV^-1) -- neutron v2 constant term
    //     vn2A     (MeV^-1) -- neutron v2 mass term
    //     vn30     (MeV^-2) -- neutron v3 constant term
    //     vn3A     (MeV^-2) -- neutron v3 mass term
    //     vp20     (MeV^-1) -- proton v2 constant term
    //     vp2A     (MeV^-1) -- proton v2 mass term
    //     vp30     (MeV^-2) -- proton v3 constant term
    //     vp3A     (MeV^-2) -- proton v3 mass term
    //     v40      (MeV^-3) -- v4 (shared n/p)
    //
    //   Imaginary volume potential (Wv):
    //     wn10     (MeV)    -- neutron w1 constant term
    //     wn1A     (MeV)    -- neutron w1 mass term
    //     wp10     (MeV)    -- proton w1 constant term
    //     wp1A     (MeV)    -- proton w1 mass term
    //     w20      (MeV)    -- w2 constant term (shared n/p)
    //     w2A      (MeV)    -- w2 mass term
    //
    //   Imaginary surface potential (Wd):
    //     d10      (MeV)    -- d1 constant term
    //     d1alpha  (MeV)    -- d1 asymmetry coefficient
    //     d20      (MeV^-1) -- d2 constant term (shared n/p)
    //     d2A      (MeV^-1) -- d2 mass term
    //     d2A2     (none)   -- d2 exponential denominator
    //     d2A3     (none)   -- d2 exponential shift
    //     d30      (MeV)    -- d3 (shared n/p)
    //
    //   Real spin-orbit potential (Vso):
    //     vSO10    (MeV)    -- vso1 constant term (shared n/p)
    //     vSO1A    (MeV)    -- vso1 mass term
    //     vSO20    (MeV^-1) -- vso2 (shared n/p)
    //
    //   Imaginary spin-orbit potential (Wso):
    //     wSO10    (MeV)    -- wso1 (shared n/p)
    //     wSO20    (MeV)    -- wso2 (shared n/p)
    //
    //   Geometry:
    //     rV0      (fm)     -- real volume radius coefficient
    //     rVA      (fm)     -- real volume radius mass term
    //     aV0      (fm)     -- real volume diffuseness constant
    //     aVA      (fm)     -- real volume diffuseness mass term
    //     rD0      (fm)     -- surface abs. radius coefficient
    //     rDA      (fm)     -- surface abs. radius mass term
    //     anD0     (fm)     -- neutron surface diffuseness constant
    //     anDA     (fm)     -- neutron surface diffuseness mass term
    //     apD0     (fm)     -- proton surface diffuseness constant
    //     apDA     (fm)     -- proton surface diffuseness mass term
    //     rSO0     (fm)     -- spin-orbit radius coefficient
    //     rSOA     (fm)     -- spin-orbit radius mass term
    //     aSO0     (fm)     -- spin-orbit diffuseness (shared n/p)
    //     rC0      (fm)     -- Coulomb radius coefficient
    //     rCA      (fm)     -- Coulomb radius mass term (1st)
    //     rCA2     (fm^2)   -- Coulomb radius mass term (2nd)
    //
    // Per-nucleus overrides may be added using the nuclear PDG code as
    // the key (e.g., "1000180400" for ⁴⁰Ar). Each per-nucleus override
    // must be a COMPLETE set of all 46 parameters; partial overrides are
    // not yet supported.
    //
    // The most convenient way to specify the "opt_mod" value is to use
    // the #include syntax to pull in one of the files from the
    // data/optical_model/ directory. For example:
    //
    //   { type: "optical_model", name: "KD-global-OMP",
    //     opt_mod: #include:"optical_model_kd_global.js" },
    //
    // The example below is equivalent to the #include form above; it uses the
    // results from the original KD global fit
    // (https://doi.org/10.1016/S0375-9474(02)01321-0), which was the default
    // optical model potential in MARLEY v1. The MARLEY v2 default is the
    // central value of the KDUQ Federal evaluation
    // (https://doi.org/10.1103/PhysRevC.107.014602), whose parameter values are
    // stored in data/optical_model/optical_model_kduq_federal_cv.js.
    { type: "optical_model", name: "KD-global-OMP",
      opt_mod: {
        Default: {
          v10: 59.30, v1A: 0.024, v1alpha: 21.0,
          vn20: 0.007228, vn2A: 1.48e-6,
          vn30: 1.994e-5, vn3A: 2.00e-8,
          vp20: 7.067e-3, vp2A: 4.23e-6,
          vp30: 1.729e-5, vp3A: 1.136e-8,
          v40: 7e-9,
          wn10: 12.195, wn1A: 0.0167,
          wp10: 14.667, wp1A: 9.629e-3,
          w20: 73.55, w2A: 0.0795,
          d10: 16.0, d1alpha: 16.0,
          d20: 0.0180, d2A: 3.802e-3,
          d2A2: 8.00, d2A3: 156.0, d30: 11.5,
          vSO10: 5.922, vSO1A: 3.00e-3, vSO20: 4.00e-3,
          wSO10: -3.10, wSO20: 160.0,
          rV0: 1.3039, rVA: 0.4054,
          aV0: 0.6778, aVA: 1.487e-4,
          rD0: 1.3424, rDA: 1.585e-2,
          anD0: 0.5446, anDA: 1.656e-4,
          apD0: 0.5187, apDA: 5.205e-4,
          rSO0: 1.1854, rSOA: 0.6470, aSO0: 0.590,
          rC0: 1.198, rCA: 0.697, rCA2: 12.994,
        },
      },
    },

    // STRENGTH VARIATION WEIGHT CALCULATOR
    //
    // A "strength_variation" weight calculator varies nuclear matrix elements
    // for discrete allowed transitions according to their experimental
    // uncertainties. Event weights are computed as the ratio of the varied
    // strength to the nominal strength for the discrete nuclear transition that
    // occurred. Only events with primary interactions corresponding to the
    // (anti-)ν CC (Discrete) and NC (Discrete) reaction types are processed;
    // events involving other reaction types are always assigned unit weight.
    //
    // The required "reaction_file" key must name one of the reaction input
    // files listed in the "reactions" array of the original generation
    // configuration. For the reweight command, that generation configuration
    // is automatically restored from the MARLEY.JSONconfig attribute saved
    // in the first input event file; this setting simply identifies which
    // reaction channel should be varied.
    //
    // Three settings for the "mode" key are available; see the annotated
    // "marley generate" configuration (examples/config/annotated.js) for
    // details beyond the information given here.

    // multisim mode (default): uncorrelated random Gaussian variations.
    // Uses a dimidiated (bifurcated) Gaussian probability density function
    // to draw random variations of each tabulated matrix element according
    // to its (possibly asymmetric) uncertainty. The required "num_variations"
    // key specifies the number of independent random draws, each producing
    // one output event weight named "<name>_<j>". An optional "seed" key
    // sets the RNG seed for the variations (default 0), independently of
    // the main MARLEY seed. The "sigma_factor" key is not allowed in this
    // mode.
    { type: "strength_variation", name: "MultisimVar",
      num_variations: 100,
      reaction_file: "ve40ArCC_Bhattacharya2009-Discrete.react",
      seed: 54321 },

    // shift mode: deterministic ±kσ shifts of all matrix elements together. The
    // required "sigma_factor" key may be a single positive number or a JSON
    // array of positive numbers. Two event weights are computed per
    // sigma_factor value, with labels of the form "<name>-up@<k>" and
    // "<name>-down@<k>". The "num_variations" and "seed" keys are not allowed
    // in this mode.
    { type: "strength_variation", name: "SysShift",
      mode: "shift",
      sigma_factor: [0.5, 1.0, 2.0],
      reaction_file: "ve40ArCC_Bhattacharya2009-Discrete.react" },

    // unisim mode: deterministic ±kσ shift of one matrix element at a time,
    // leaving all others at their nominal values. Two event weights are
    // computed per matrix element per sigma_factor value, with labels of the
    // form "<name>-me<m>_up@<k>" and "<name>-me<m>_down@<k>". Here <m> is the
    // zero-based index of the varied matrix element. The "num_variations" and
    // "seed" keys are not allowed in this mode.
    { type: "strength_variation", name: "UniShift",
      mode: "unisim",
      sigma_factor: 1.0,
      reaction_file: "ve40ArCC_Bhattacharya2009-Discrete.react" },
  ],

  // REWEIGHT COMMAND SETTINGS (optional)
  //
  // The "reweight" JSON object controls execution of the "marley reweight"
  // command. Its "output" JSON array uses the same format as the
  // "generate.output" array described in examples/config/annotated.js,
  // with two differences: the default output file name is
  // "reweighted_events.hepmc3", and only "overwrite" mode is supported
  // (the "resume" mode is not available for the reweight command).
  reweight: {

    output: [ { file: "reweighted_events.hepmc3", format: "ascii",
      mode: "overwrite" } ],
  },

} // A closing curly brace should appear at the end of the file
