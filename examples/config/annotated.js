// Example MARLEY job configuration file
// Steven Gardiner <gardiner@fnal.gov>
// Revised 11 July 2026 for MARLEY 2.0.0
//
// INTRODUCTION
//
// The MARLEY command-line executable is configured using a JSON-like file
// format. While the file format is quite similar to standard JSON (see
// http://www.json.org/ for a full description), MARLEY job configuration files
// differ from JSON files in the following ways:
//
//   - Single-word keys (no whitespace) may be given without surrounding
//     double quotes
//
//   - C++-style comments // and /* */ are allowed anywhere in the file.
//
//   - A trailing comma is allowed at the end of JSON objects and arrays
//
//  The JSON parser included with MARLEY will print an error message
//  if parsing the job configuration file fails. This message will
//  provide some (hopefully useful) guidance in troubleshooting
//  formatting mistakes. To simplify writing their own configuration files,
//  users are encouraged to copy the examples (particularly
//  "examples/COPY_ME.js") and modify them to suit their needs.
//
//  A file extension of ".js" is recommended for MARLEY configuration files
//  because typical syntax highlighting settings for JavaScript work well with
//  the configuration file format.
//
{ // An opening curly brace begins the configuration file content

  // RANDOM NUMBER SEED (optional)
  //
  // The "seed" key provides a 64-bit unsigned integer (between 0
  // and 2^64 - 1, inclusive) that will be used to seed MARLEY's random
  // number generator.
  //
  // If this key is omitted, MARLEY will use the system time since the
  // Unix epoch as its random number seed.
  seed: 123456,

  // INCIDENT NEUTRINO DIRECTION (optional)
  //
  // The "direction" JSON object stores a 3-vector that represents the
  // direction of the incident neutrinos. Note that this vector does not need
  // to be normalized to unity, but at least one element *must* be nonzero.
  //
  // The "x", "y", and "z" keys give the Cartesian components of the direction
  // 3-vector. If the "direction" JSON object is omitted, then x = 0.0, y =
  // 0.0, z = 1.0 will be assumed.
  //
  // An isotropic neutrino direction may be randomly sampled for each event by
  // using the following configuration:
  //
  // direction: "isotropic",
  //
  // In this example, incident neutrinos travel in the +z direction.
  direction: { x: 0.0, y: 0.0, z: 1.0 },

  // TARGET SPECIFICATION (optional)
  //
  // The nuclidic composition of the material illuminated by the incident
  // neutrinos in a MARLEY simulation may be specified using the "target"
  // JSON object. This object defines two arrays, which must have equal sizes.
  // The "nuclides" array contains one or more nuclear PDG codes with one entry
  // per distinct nuclide present in the target material. The "atom_fractions"
  // array, as its name suggests, contains the corresponding atom fractions.
  // Elements of the "atom_fractions" array will automatically be normalized
  // to sum to unity if this is not done already by the user. Negative
  // elements in this array will trigger an error message from MARLEY.
  //
  // If the "target" object is omitted from the job configuration file,
  // then every nuclide that appears in the initial state of at least one
  // configured reaction (see the REACTION INPUT FILES(S) section below)
  // will be assumed to be present with equal abundance.
  target: {
    nuclides: [ 1000180400 ],
    atom_fractions: [ 1.0 ],
  },

  // REACTION INPUT FILE(S) (required)
  //
  // For simulating neutrino-nucleus reactions, MARLEY relies on tables of
  // precomputed nuclear matrix elements to compute neutrino reaction cross
  // sections. The nuclear matrix elements to use in any particular simulation
  // job are specified using a JSON array stored under the "reactions" key.
  // Each entry in the array is the name of a data file that contains a set of
  // nuclear matrix elements to use when calculating cross sections for a
  // particular scattering mode on a particular target nuclide. Note that the
  // file names must be simple strings; MARLEY configuration files currently do
  // not support the use of environment variables, bash globs, etc.

  // MARLEY also uses an input file to configure neutrino-electron
  // elastic scattering (ES) reactions. In this case, the file provides a list
  // of the nuclides for which the ES process should be simulated.

  // Although multiple reaction input files may be available for a particular
  // process, only a single configuration is allowed for any particular
  // combination of target nucleus and scattering mode. Conflicting
  // configurations will be ignored in favor of the one that appears in the
  // first relevant reaction input file listed in the "reactions" array.
  //
  // In this version of MARLEY, the main process available to be simulated is
  // charged-current scattering of electron neutrinos on 40Ar. Reaction input
  // files are stored in the folder data/react/.
  //
  // Two reaction input files bundled with MARLEY for processes other than
  // CC neutrino-40Ar scattering are:
  //
  //   - ES.react: Enables simulation of neutrino-electron elastic scattering
  //               on a 40Ar atomic target
  //
  //   - CEvNS40Ar.react: Enables simulation of coherent elastic
  //                      neutrino-nucleus scattering on a 40Ar target. The
  //                      q^2 dependence of the nuclear form factor is
  //                      neglected.
  //
  // The full path to the file does not need to be given in each element of the
  // "reactions" JSON array. After searching in the working directory, MARLEY's
  // default behavior is to search for reaction input files in the directories
  // ${MARLEY}/data, ${MARLEY}/data/react/, and ${MARLEY}/data/structure/,
  // where ${MARLEY} is the value of the MARLEY environment variable (typically
  // set by sourcing the setup_marley.sh bash script). The list of search
  // directories beyond the working directory can be changed from the default
  // by setting the MARLEY_SEARCH_PATH environment variable to a ':'-separated
  // list of directories.
  //
  // Unless a particular reaction channel is represented by a data file given
  // in the "reactions" JSON array, it will not be included in the MARLEY
  // simulation.
  //
  reactions: [ "ve40ArCC_Bhattacharya2009.react", "ES.react" ],

  // ---- ADVANCED OPTIONS (all optional) ----
  //
  // The keys described in this section are advanced configuration options
  // that most users will not need to adjust. Sensible defaults are used
  // whenever these keys are omitted from the job configuration file.

  // COULOMB CORRECTION METHOD (optional)
  //
  // The "coulomb_mode" key selects the method used to compute Coulomb
  // corrections for charged-current nuclear reactions. Valid values are:
  //
  //   - "none": No Coulomb correction applied
  //   - "Fermi": Use the Fermi function
  //   - "EMA": Use the effective momentum approximation (EMA)
  //   - "MEMA": Use a modified version of the EMA
  //   - "Fermi-EMA": Interpolate between the Fermi function and the EMA
  //   - "Fermi-MEMA": Interpolate between the Fermi function and the MEMA
  //                   (this is the default)
  //
  // coulomb_mode: "Fermi-MEMA",

  // NUCLEON AND NUCLEAR FORM FACTORS (optional)
  //
  // The "form_factors" key controls the parameterizations used for nucleon
  // (Sachs and axial) and nuclear form factors. It may be set to an object:
  //
  //   form_factors: {
  //     sachs_model: "bbba05",    // "trivial", "dipole", or "bbba05"
  //     axial_model: "dipole",    // "trivial" or "dipole"
  //     nuclear_model: "klein",   // "trivial", "helm", or "klein"
  //
  //     // When nuclear_model is "klein", an optional sub-configuration
  //     // may be provided:
  //     // nucl_options: { adapted: false },
  //   },
  //
  // The default configuration shown above will be used when the
  // "form_factors" key is omitted.
  //
  // Alternatively, this key may be set to the string "allowed", "AA", or
  // "aa" to use the allowed approximation (trivial form factors for all
  // three categories).

  // NUCLEAR DE-EXCITATIONS (optional)
  //
  // Use a boolean value to enable or disable simulation of nuclear
  // de-excitations for all reactions. The default is true.
  //
  // do_deexcitations: true,

  // SUB-CONTINUUM MODE (optional)
  //
  // The "sub_continuum_mode" key controls how cross-section strength
  // that falls below the unbound threshold is handled. Valid values are:
  //
  //   - "ignore": Cross-section strength below threshold is discarded
  //   - "mirror": Strength is mirrored from above the threshold
  //   - "accumulate": Strength accumulates at the threshold (this is
  //                   the default)
  //
  // sub_continuum_mode: "accumulate",

  // OPTICAL MODEL PARAMETERS (optional)
  //
  // The "opt_mod" key provides a custom configuration of nuclear optical
  // model parameters, overriding the defaults. The value should be a JSON
  // object whose format matches the optical model configuration used by
  // the MARLEY structure database.
  //
  // opt_mod: { ... },

  // ANGULAR MOMENTUM CUTOFFS (optional)
  //
  // The "fragment_lmax" key sets the maximum orbital angular momentum
  // quantum number to consider when computing fragment decay widths
  // (default: 2). The "gamma_lmax" key sets the maximum multipolarity
  // for gamma-ray decay widths (default: 2). Both values are integers;
  // gamma_lmax must be >= 1.
  //
  // fragment_lmax: 2,
  // gamma_lmax: 2,

  // ENERGY PDF MAXIMUM (optional)
  //
  // If MARLEY has difficulty automatically finding the maximum of the
  // neutrino energy probability density function (this can happen for
  // unusual user-defined spectra), you may provide your own estimate via
  // the "energy_pdf_max" key. The value should be an energy in MeV.
  //
  // energy_pdf_max: 50.0,

  // NEUTRINO SOURCE SPECIFICATION (required)
  //
  // The "source" JSON object describes the incident neutrino spectrum. The
  // input spectrum should not be weighted by a reaction cross section because
  // MARLEY will weight the incident spectrum using its own cross section
  // model during the simulation. The flux units used to specify user-defined
  // spectra will be ignored by MARLEY (the source spectrum will be
  // renormalized to unity internally to produce a probability density
  // function), but all neutrino energies should be given in MeV.
  //
  // The "type" key owned by the source specification describes the
  // incident spectrum and determines the other parameters that must be
  // specified. The currently allowed values for the "type" key are
  // given in the following table (with synonyms separated by commas)
  //
  //   Spectrum description                   Allowed "type" values
  //   --------------------                   ---------------------
  //
  //   Fermi-Dirac                            "fd", "fermi-dirac",
  //                                          "fermi_dirac"
  //
  //   "Alpha-fit"                            "af", "alpha", "alpha-fit"
  //   (see, e.g.,
  //   https://arxiv.org/abs/astro-ph/0208035)
  //
  //   "Beta-fit"                             "bf", "beta", "beta-fit"
  //   (see, e.g.,
  //   https://arxiv.org/abs/1511.00806)
  //
  //   Monoenergetic                          "mono", "monoenergetic"
  //
  //   Muon decay-at-rest                     "dar", "decay-at-rest"
  //   (ve and vu only)
  //
  //   User-defined histogram                 "hist", "histogram"
  //
  //   User-defined probability               "grid"
  //   density function evaluated via
  //   interpolation on a set of grid
  //   points
  //
  //   ROOT TH1                               "th1"
  //
  //   ROOT TGraph                            "tgraph"
  //
  // All of the source types also require the use of the "neutrino" key,
  // which specifies the species of neutrino emitted by the source. Valid
  // values for the "neutrino" key are neutrino PDG codes (±12, ±14, ±16)
  // and the strings "ve", "vebar", "vu", "vubar", "vt", and "vtbar".
  //
  // Although MARLEY's algorithm for sampling from arbitrary spectra is
  // reasonably robust, user-defined spectra with widely-separated tight
  // peaks or other unusual shapes may experience problems. Users should
  // perform some simple verification simulations when using a "hist",
  // "grid", "th1", or "tgraph" neutrino source. For assistance with
  // this testing or to report bugs, please contact the MARLEY developers
  // (support@marleygen.org).
  //
  // Examples of each of the allowed source specifications are shown below
  //
  // FERMI-DIRAC
  //
  //  source: {
  //    type: "fermi-dirac",
  //    neutrino: "ve",
  //    Emin: 0,           // Minimum neutrino energy (MeV)
  //    Emax: 60,          // Maximum neutrino energy (MeV)
  //    temperature: 3.5,  // Temperature (MeV)
  //    eta: 4             // Pinching parameter (dimensionless, default 0)
  //  },
  //
  // "ALPHA FIT"
  //
  //  source: {
  //    type: "alpha-fit",
  //    neutrino: "ve",
  //    Emin: 0,           // Minimum neutrino energy (MeV)
  //    Emax: 60,          // Maximum neutrino energy (MeV)
  //    Emean: 15,         // Mean neutrino energy (MeV)
  //    alpha: 2.0,        // Pinching parameter (dimensionless, default 2.0)
  //  },
  //
  // "BETA FIT"
  //
  //  source: {
  //    type: "beta-fit",
  //    neutrino: "ve",
  //    Emin: 0,           // Minimum neutrino energy (MeV)
  //    Emax: 60,          // Maximum neutrino energy (MeV)
  //    Emean: 15,         // Mean neutrino energy (MeV)
  //    beta: 3.0,         // Pinching parameter (dimensionless, default 4.5)
  //  },
  //
  //  MONOENERGETIC
  //
  //  source: {
  //    type: "monoenergetic",
  //    neutrino: "ve",
  //    energy: 10,        // Neutrino energy (MeV)
  //  },
  //
  //  MUON DECAY-AT-REST
  //
  //  source: {
  //    type: "decay-at-rest",
  //    neutrino: "ve",
  //  },
  //
  //  HISTOGRAM
  //
  //  source: {
  //    type: "histogram",
  //    neutrino: "ve",
  //    E_bin_lefts: [ 10., 20., 30. ],   // Low edges of energy bins (MeV)
  //    weights: [ 0.2, 0.5, 0.3 ],       // Bin weights (dimensionless)
  //    Emax: 40.,                        // Upper edge of the final bin (MeV)
  //  },
  //
  //  Within a histogram bin, energies are sampled uniformly on the
  //  half-open interval [ Ebin_left, Ebin_right ).
  //
  //  GRID
  //
  //  source: {
  //    type: "grid",
  //    neutrino: "ve",
  //    energies: [ 10., 15., 20. ],   // Energy grid points (MeV)
  //
  //    prob_densities: [ 0., 1., 0. ],  // Probability densities
  //                                     // (dimensionless, do not need to be
  //                                     //  normalized to unity by the user)
  //
  //    rule: "linlin",                // Interpolation rule ("linlin" default)
  //  },
  //
  //  The allowed values of the "rule" key are "linlin" (linear-linear
  //  interpolation), "loglog" (log-log interpolation), "linlog" (linear
  //  in energy, logarithmic in probability density), and "loglin"
  //  (logarithmic in energy, linear in probability density).
  //
  //
  //  TH1
  //
  //  source: {
  //    type: "th1",
  //    neutrino: "ve",
  //    tfile: "my_root_file.root",  // Name of the ROOT file containing
  //                                 // the TH1 object
  //
  //    namecycle: "MyHist",         // Name under which the TH1 object
  //                                 // appears in the file (used to
  //                                 // retrieve it)
  //  },
  //
  //  TGRAPH
  //
  //  source: {
  //    type: "tgraph",
  //    neutrino: "ve",
  //    tfile: "my_root_file.root",  // Name of the ROOT file containing
  //                                 // the TGraph object
  //
  //    namecycle: "MyGraph",        // Name of the TGraph object (used to
  //                                 // retrieve it from the ROOT file)
  //  },
  //
  //
  // In this example configuration file, we've chosen a monoenergetic source.
  //
  source: {
    neutrino: "ve",        // The source produces electron neutrinos
    type: "monoenergetic",
    energy: 15.0,          // MeV

    // WEIGHT FLUX (optional)
    //
    // By default, MARLEY weights the incident neutrino spectrum by the
    // total reaction cross section(s) when generating events. Set this
    // boolean key to false to sample neutrino energies directly from the
    // unweighted source spectrum instead.
    //
    // weight_flux: true,
  },

  // EVENT WEIGHT CALCULATORS (optional)
  //
  // The "weights" JSON array configures one or more weight calculators
  // that assign additional event weights during generation. Each entry is
  // a JSON object with at least a "type" key and a "name" key.
  //
  // Available weight calculator types:
  //
  //   - "trivial": Returns unit weight for every event.
  //                Required keys: type, name
  //
  //   - "optical_model": Recomputes Hauser-Feshbach decay widths using
  //                      an alternative set of optical model parameters.
  //                      Required keys: type, name, opt_mod
  //
   //   - "strength_variation": Varies nuclear matrix element strengths
   //                            according to their experimental
   //                            uncertainties. Two modes are available:
   //
   //     "multisim" (default): Uses a dimidiated (bifurcated) Gaussian
   //                            probability density function to draw
   //                            random variations. Event weights are
   //                            computed as the ratio of the varied
   //                            strength to the nominal strength for
   //                            events corresponding to discrete nuclear
   //                            transitions.
   //                            Required keys: type, name, num_variations,
   //                            reaction_file
   //                            Optional keys: mode, seed (default 0;
   //                            separate from the main MARLEY RNG seed set
   //                            by the top-level "seed" key)
   //
   //     "sigma_shift":       Applies a deterministic +k or -k sigma shift
   //                            to each matrix element strength, where k
   //                            is given by sigma_factor and sigma is the
   //                            experimental uncertainty. Intended for
   //                            systematic variation studies. Two weight
   //                            calculators are created per sigma_factor
   //                            value: "<name>-up@<k>" and
   //                            "<name>-down@<k>".
   //                            Required keys: type, name, reaction_file,
   //                            sigma_factor (scalar or array of positive
   //                            numbers)
   //                            Optional keys: mode
  //
   // weights: [
   //   { type: "trivial", name: "MyWeight" },
   //   { type: "optical_model", name: "OMP", opt_mod: { ... } },
   //
   //   // Multisim mode (default): random Gaussian variations
   //   { type: "strength_variation", name: "MultisimVar",
   //     num_variations: 100,
   //     reaction_file: "my_reaction.react", seed: 12345 },
   //
   //   // Sigma_shift mode: deterministic systematic shifts
   //   { type: "strength_variation", name: "SysShift",
   //     mode: "sigma_shift",
   //     sigma_factor: [0.5, 1.0, 2.0],
   //     reaction_file: "my_reaction.react" },
   // ],

  // EXECUTABLE SETTINGS (optional)
  //
  // The entries within the executable_settings JSON object are used to
  // control the marley command-line executable. They are ignored if
  // the job configuration file is used to initialize MARLEY outside of that
  // context (e.g., within a Geant4 application that links to the MARLEY
  // shared libraries).
  executable_settings: {

    // EVENT COUNT (optional)
    //
    // Specifies the number of events to produce before terminating the
    // program. The JSON parser expects this entry to be an integer literal,
    // so scientific notation is not currently allowed.
    //
    // If this key is omitted, a value of 1000 will be assumed.
    events: 100000,

    // STATUS UPDATE INTERVAL (optional)
    //
    // Specifies the number of events between status display updates.
    // The value must be a positive integer. The default is 100.
    //
    // status_update_interval: 100,

    // EVENT OUTPUT (optional)
    //
    // The "output" JSON array contains a list of JSON objects representing
    // zero or more output streams that will receive the events generated
    // by the marley command-line executable
    //
    // Each entry is a JSON object with the following keys:
    //
    //   - file: The name of a file that will store the generated events.
    //           Streaming events to stdout or stderr is not currently
    //           supported.
    //
    //   - format: The format to use when storing the events in the file.
    //             Valid values are "ascii" and "root".
    //             Details about the format options are given below.
    //
    //   - mode: The file I/O mode to use when writing to this file. Valid
    //           values are "overwrite" (erase any previously existing file
    //           contents) and "resume". If the "resume" mode
    //           is chosen, the generator will restore its previous state from
    //           an incomplete run (e.g., a run that was interrupted by the
    //           user via ctrl+C) that was saved to the output file and
    //           continue from where it left off.
    //
    //           At most one output file may use the "resume" mode. If
    //           multiple output file entries specify mode: "resume",
    //           then the marley executable will report an error.
    //
    //   - force: Boolean value used only for the "overwrite" mode.
    //            If the output file already exists and force is not
    //            set to true, the marley executable will prompt the
    //            user before overwriting it. If the user declines,
    //            execution stops with an error message that suggests
    //            how to adjust the configuration (set "force" to
    //            true to overwrite automatically, or use
    //            "mode": "resume" to append to the existing file).
    //            If this key is omitted, a value of false is assumed.
    //
    // The allowed output file formats are
    //
    //   - "ascii": Events are stored in the standard HepMC3 text format
    //              with full floating-point precision. This format uses the
    //              HepMC3::WriterAscii class for output and the
    //              HepMC3::ReaderAscii class for input. The output is
    //              compliant with the NuHepMC v1.0.0 standard for neutrino
    //              event generators, making it readable by any
    //              HepMC3-compatible application.
    //
    //   - "root": Events are stored in a compressed ROOT TTree named
    //             "MARLEY_event_tree" using the HepMC3 GenEventData POD
    //             structure. Because GenEventData is a plain-old-data struct
    //             with an associated ROOT dictionary distributed with HepMC3,
    //             the output file can be read without loading MARLEY-specific
    //             shared libraries. This format is only available if MARLEY
    //             has been built with ROOT support.
    //
  // If this key is omitted, then the following configuration
  // is assumed:
  //
  // output: [ { file: "events.hepmc3", format: "ascii", mode: "overwrite",
  //             force: false } ]
  //
  output: [ { file: "events.hepmc3", format: "ascii", mode: "overwrite" } ],
  },

  // LOGGER CONFIGURATION (separate file)
  //
  // MARLEY's diagnostic message output is configured separately from the
  // job configuration file via the file data/config/logger.js (located
  // using the MARLEY environment variable). This file uses a JSON format
  // to configure output streams, severity levels, and logging categories.
  // See the comments in that file for details.

} // A closing curly brace should appear at the end of the file
