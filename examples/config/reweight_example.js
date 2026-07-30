// Example MARLEY job configuration file for the "marley reweight" command
// Steven Gardiner <gardiner@fnal.gov>
// Revised 19 July 2026 for MARLEY 2.0.0
//
// INTRODUCTION
//
// The "marley reweight" command assigns systematic event weights to a
// previously generated HepMC3 event file. It reads a configuration file
// (see below) that specifies one or more weight calculators, applies them
// to each event in the input file, and writes the reweighted events to an
// output file.
//
// The configuration format for the "marley reweight" command uses the same
// JSON-like syntax as the other MARLEY commands. The "weights" top-level
// key specifies the weight calculator(s) to apply, while the optional
// "reweight" section may be used to configure output settings. All other
// parameters needed by the Generator are automatically restored from the
// "MARLEY.JSONconfig" attribute saved in the input event file.
//
{ // An opening curly brace begins the configuration file content

  // EVENT WEIGHT CALCULATORS (required)
  //
  // The "weights" JSON array configures one or more weight calculators
  // that assign additional event weights. Each entry is a JSON object
  // with at least a "type" key and a "name" key.
  //
  // The following weight calculator types are available:
  //
  //   "trivial"            Returns unit weight for every event
  //   "optical_model"      Recomputes Hauser-Feshbach decay widths using
  //                        an alternative optical model
  //   "strength_variation" Varies nuclear matrix element strengths
  //                        according to their experimental uncertainties
  //
  weights: [

    //////// TRIVIAL ////////
    //
    // A trivial weight calculator assigns a weight of exactly 1.0 to
    // every event. Its main use is to label a weight column in the
    // output event file with a descriptive name so that downstream
    // analysis code can refer to it.
    //
    // Parameters:
    //   type  (string, required)  Must be "trivial"
    //   name  (string, required)  Name used to label the weight column
    //
    { type: "trivial", name: "MyTrivialWeight" },

    //////// OPTICAL MODEL ////////
    //
    // An optical model weight calculator reconstructs the Hauser-Feshbach
    // decay widths using a user-specified alternative set of optical model
    // parameters and computes the weight as the ratio of the alternative
    // width to the original simulation width.
    //
    // Parameters:
    //   type    (string, required)  Must be "optical_model"
    //   name    (string, required)  Name used to label the weight column
    //   opt_mod (object, required)  Custom optical model parameters
    //
    // The "opt_mod" object contains a flat map of string keys to JSON
    // sub-objects (see details below). At minimum, the key "Default"
    // must be present with 46 double-valued parameters defining the
    // Koning-Delaroche optical model potential (Koning & Delaroche,
    // Nucl. Phys. A 713, 231 (2003), https://doi.org/10.1016/S0375-9474(02)01321-0).
    // Additional per-nucleus overrides may be specified using the PDG code
    // of the target nuclide as the key (e.g., "1000180400" for 40Ar).
    //
    // The 46 required parameters for each sub-object (shown here with
    // illustrative values; see data/optical_model/optical_model_kduq_federal_cv.js
    // for the actual default MARLEY v2.0.0 values) are:
    //
    //   Real volume potential (Vv):
    //     v10      (MeV)     55.0        -- v1 constant term
    //     v1A      (MeV)      0.02       -- v1 mass-number coefficient
    //     v1alpha  (MeV)     11.0        -- v1 asymmetry coefficient
    //     vn20     (MeV^-1)   0.005      -- neutron v2 constant term
    //     vn2A     (MeV^-1)   0.000001   -- neutron v2 mass term
    //     vn30     (MeV^-2)   0.00001    -- neutron v3 constant term
    //     vn3A     (MeV^-2)   0.00000001 -- neutron v3 mass term
    //     vp20     (MeV^-1)   0.006      -- proton v2 constant term
    //     vp2A     (MeV^-1)   0.000004   -- proton v2 mass term
    //     vp30     (MeV^-2)   0.00001    -- proton v3 constant term
    //     vp3A     (MeV^-2)   0.00000003 -- proton v3 mass term
    //     v40      (MeV^-3)  -0.00000001 -- v4 (shared n/p)
    //
    //   Imaginary volume potential (Wv):
    //     wn10     (MeV)     20.0        -- neutron w1 constant term
    //     wn1A     (MeV)      0.01       -- neutron w1 mass term
    //     wp10     (MeV)     18.0        -- proton w1 constant term
    //     wp1A     (MeV)      0.03       -- proton w1 mass term
    //     w20      (MeV)    100.0        -- w2 constant term (shared n/p)
    //     w2A      (MeV)      0.05       -- w2 mass term
    //
    //   Imaginary surface potential (Wd):
    //     d10      (MeV)     20.0        -- d1 constant term
    //     d1alpha  (MeV)     12.0        -- d1 asymmetry coefficient
    //     d20      (MeV^-1)   0.03       -- d2 constant term (shared n/p)
    //     d2A      (MeV^-1)   0.003      -- d2 mass term
    //     d2A2     (none)     8.0        -- d2 exponential denominator
    //     d2A3     (none)   250.0        -- d2 exponential shift
    //     d30      (MeV)     16.0        -- d3 (shared n/p)
    //
    //   Real spin-orbit potential (Vso):
    //     vSO10    (MeV)      5.5        -- vso1 constant term (shared n/p)
    //     vSO1A    (MeV)     -0.0005     -- vso1 mass term
    //     vSO20    (MeV^-1)   0.005      -- vso2 (shared n/p)
    //
    //   Imaginary spin-orbit potential (Wso):
    //     wSO10    (MeV)     -3.5        -- wso1 (shared n/p)
    //     wSO20    (MeV)    250.0        -- wso2 (shared n/p)
    //
    //   Geometry:
    //     rV0      (fm)      1.29        -- real volume radius coefficient
    //     rVA      (fm)      0.41        -- real volume radius mass term
    //     aV0      (fm)      0.67        -- real volume diffuseness constant
    //     aVA      (fm)     -0.0002      -- real volume diffuseness mass term
    //     rD0      (fm)      1.37        -- surface abs. radius coefficient
    //     rDA      (fm)      0.02        -- surface abs. radius mass term
    //     anD0     (fm)      0.56        -- neutron surface diffuseness constant
    //     anDA     (fm)     -0.00006     -- neutron surface diffuseness mass term
    //     apD0     (fm)      0.49        -- proton surface diffuseness constant
    //     apDA     (fm)      0.001       -- proton surface diffuseness mass term
    //     rSO0     (fm)      1.21        -- spin-orbit radius coefficient
    //     rSOA     (fm)      0.66        -- spin-orbit radius mass term
    //     aSO0     (fm)      0.57        -- spin-orbit diffuseness (shared n/p)
    //     rC0      (fm)      1.21        -- Coulomb radius coefficient
    //     rCA      (fm)      0.77        -- Coulomb radius mass term (1st)
    //     rCA2     (fm^2)   13.0         -- Coulomb radius mass term (2nd)
    //
    //   Numerov integration:
    //     step_size (fm)     0.1         -- Numerov step size (optional)
    //
    // The example below uses modified values for illustration.
    //
    // Per-nucleus overrides may be added using the PDG code of the
    // nuclide as the key (e.g., "1000180400" for 40Ar). Each per-nucleus
    // override must be a COMPLETE set of all parameters; partial overrides
    // are not supported.
    //
    { type: "optical_model", name: "OMP_Alt",
      opt_mod: {
        Default: {
          v10: 55.0, v1A: 0.02, v1alpha: 11.0,
          vn20: 0.005, vn2A: 0.000001,
          vn30: 0.00001, vn3A: 0.00000001,
          vp20: 0.006, vp2A: 0.000004,
          vp30: 0.00001, vp3A: 0.00000003,
          v40: -0.00000001,
          wn10: 20.0, wn1A: 0.01,
          wp10: 18.0, wp1A: 0.03,
          w20: 100.0, w2A: 0.05,
          d10: 20.0, d1alpha: 12.0,
          d20: 0.03, d2A: 0.003,
          d2A2: 8.0, d2A3: 250.0, d30: 16.0,
          vSO10: 5.5, vSO1A: -0.0005, vSO20: 0.005,
          wSO10: -3.5, wSO20: 250.0,
          rV0: 1.29, rVA: 0.41,
          aV0: 0.67, aVA: -0.0002,
          rD0: 1.37, rDA: 0.02,
          anD0: 0.56, anDA: -0.00006,
          apD0: 0.49, apDA: 0.001,
          rSO0: 1.21, rSOA: 0.66, aSO0: 0.57,
          rC0: 1.21, rCA: 0.77, rCA2: 13.0,
          step_size: 0.1,
        },
      },
    },

    //////// STRENGTH VARIATION ////////
    //
    // A strength variation weight calculator varies nuclear matrix element
    // strengths according to their experimental uncertainties (as encoded
    // in the reaction input file). Three operation modes are available:
    // "multisim", "shift", and "unisim".
    //
    // Common parameters for all modes:
    //   type           (string, required)  Must be "strength_variation"
    //   name           (string, required)  Base name for weight columns
    //   reaction_file  (string, required)  Path to a reaction input file
    //                     that contains discrete nuclear transitions with
    //                     experimental strength uncertainties
    //
    // Multisim mode (default):
    //   Draws strength variations from a dimidiated (bifurcated) Gaussian
    //   PDF. Each random draw creates one weight calculator instance named
    //   "<name>_<idx>". Several instances are typically created to build
    //   a distribution of weights.
    //
    //   Parameters:
    //     mode            (string, optional)  "multisim" (default)
    //     num_variations  (int, required)     Number of random variations
    //     seed            (int, optional)     RNG seed (default 0; separate
    //                       from the main MARLEY seed)
    //
    //   Note: sigma_factor is NOT allowed in multisim mode.
    //
    { type: "strength_variation", name: "MultisimVar",
      num_variations: 100,
      reaction_file: "ve40ArCC_Bhattacharya2009-Discrete.react",
      seed: 54321 },

    // Sigma_shift mode:
    //   Applies a deterministic +k or -k sigma shift to every matrix
    //   element simultaneously, where k = sigma_factor. Two calculators
    //   are created per sigma_factor value: "<name>-up@<k>" and
    //   "<name>-down@<k>". Intended for systematic variation studies.
    //
    //   Parameters:
    //     mode          (string, required)   "shift"
    //     sigma_factor  (scalar or array,    Positive number(s)
    //                    required)           defining the shift size
    //
    //   Note: num_variations and seed are NOT allowed in shift mode.
    //
    { type: "strength_variation", name: "SysShift",
      mode: "shift",
      sigma_factor: [0.5, 1.0, 2.0],
      reaction_file: "ve40ArCC_Bhattacharya2009-Discrete.react" },

    // Unisim mode:
    //   Applies a deterministic +k or -k sigma shift to one individual
    //   matrix element at a time, leaving all others at their nominal
    //   values. Two calculators are created per sigma_factor value per
    //   matrix element: "<name>-me<m>_up@<k>" and
    //   "<name>-me<m>_down@<k>", where <m> is the zero-based index of
    //   the varied matrix element. Intended for univariate sensitivity
    //   studies to identify which transitions drive the overall uncertainty.
    //
    //   Parameters:
    //     mode          (string, required)   "unisim"
    //     sigma_factor  (scalar or array,    Positive number(s)
    //                    required)           defining the shift size
    //
    //   Note: num_variations and seed are NOT allowed in unisim mode.
    //
    { type: "strength_variation", name: "UniShift",
      mode: "unisim",
      sigma_factor: 1.0,
      reaction_file: "ve40ArCC_Bhattacharya2009-Discrete.react" },
  ],

  // REWEIGHT COMMAND SETTINGS (optional)
  //
  // Output settings for the "marley reweight" command are contained within
  // the "reweight" JSON object. This follows the same pattern as the
  // "generate" section used by the "marley generate" command.
  reweight: {

    // OUTPUT CONFIGURATION (optional)
    //
    // The "output" JSON array follows the same format as the one used by
    // the "marley generate" command's "generate.output" key.
    // Each entry is a JSON object with the following keys:
    //
    //   file    (string, optional)  Output file name. Default for reweight:
    //                               "reweighted_events.hepmc3"
    //
    //   format  (string, optional)  Output format. Valid values are "ascii"
    //                               and "root". Default: "ascii"
    //
    //   mode    (string, optional)  File I/O mode. Only "overwrite" is
    //                               allowed for the reweight command.
    //                               Default: "overwrite"
    //
    //   force   (bool,   optional)  Overwrite without prompting. If false
    //                               and the output file exists, the program
    //                               will ask before overwriting.
    //                               Default: false
    //
    // The "resume" mode is not supported by the reweight command. When
    // multiple output entries are specified, each must use mode "overwrite".
    //
    output: [ { file: "reweighted_events.hepmc3", format: "ascii",
      mode: "overwrite" } ],
  },

} // A closing curly brace should appear at the end of the file
