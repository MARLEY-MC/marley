// Example MARLEY configuration file
// S. Gardiner <sjgardiner@ucdavis.edu>
// Revised 22 September 2018 for version 1.1.0
//
// INTRODUCTION
//
// Starting with version 1.0.0, the MARLEY command-line executable
// is configured using a JSON-like file format. While the file
// format is quite similar to standard JSON (see http://www.json.org/
// for a full description), MARLEY configuration files differ from
// JSON files in the following ways:
//
//   - Single-word keys (no whitespace) may be given without surrounding
//     double quotes
//
//   - C++-style comments // and /* */ are allowed anywhere in the file.
//
//   - A trailing comma is allowed at the end of JSON objects and arrays
//
//  The JSON parser included with MARLEY attempts to be as permissive as
//  possible. It will usually accept input even if there are some minor
//  formatting mistakes (e.g., a missing opening curly brace).
//
//  A file extension of ".js" is recommended for MARLEY configuration files
//  because typical syntax highlighting settings for Javascript work well with
//  the configuration file format.
//
{  // An opening curly brace begins the configuration file content

  // RANDOM NUMBER SEED (optional)
  //
  // The "seed" key provides a 64-bit unsigned integer (between 0
  // and 2^64 - 1, inclusive) that will be used to seed MARLEY's random
  // number generator.
  //
  // If this key is omitted, MARLEY will use the system time since the
  // Unix epoch as its random number seed.
  seed: 123456,

  // NUCLEAR STRUCTURE DATA (optional, *strongly* recommended)
  //
  // The generator uses discrete nuclear level data and gamma-ray branching
  // ratios to simulate nuclear de-excitations. The nuclear structure data
  // included with MARLEY were taken with attribution from the TALYS nuclear
  // reaction code (see structure/README.md for more information) and
  // reformatted.
  //
  // The "structure" JSON array contains a list of nuclear structure data files
  // that should be used during the simulation. Paths must be included unless
  // the data files are stored in the current working directory. Relative paths
  // are assumed to be specified relative to the current working directory.
  // MARLEY configuration files currently do not support the use of environment
  // variables, bash globs, etc.
  //
  // MARLEY will load and use data for all nuclides present in each file given
  // in this list. If MARLEY creates a final-state nucleus
  // for which no nuclear structure data are available, the code will
  // use simple gamma-ray models to simulate bound-state de-excitations.
  //
  // The three files shown in this example contain all of the structure data
  // currently recommended for use with MARLEY.
  structure: [ "../structure/z019",
               "../structure/z018",
               "../structure/z017", ],

  // REACTION DATA (required)
  //
  // The generator uses tables of nuclear matrix elements to compute neutrino
  // reaction cross sections. Each entry in the "reactions" JSON array is the
  // name of a data file that contains a set of nuclear matrix elements to use
  // when calculating cross sections for a particular reaction channel (e.g.,
  // charged current ve, charged current anti-ve, neutral current, elastic
  // scattering on atomic electrons). Although multiple nuclear matrix element
  // tables may be available for a particular channel (corresponding to
  // different data evaluations) the user should ensure that only one
  // reaction data file for each desired reaction channel is used.
  //
  // In MARLEY v1.1.1, the only available channel is charged current ve,
  // and there are three available tables of evaluated nuclear matrix elements
  // for this channel:
  //
  //   - react/ve40ArCC_Bhattacharya2009.react
  //   - react/ve40ArCC_Bhattacharya1998.react
  //   - react/ve40ArCC_Liu1998.react
  //
  // See each of these files for details about their respective nuclear
  // matrix element evaluations.
  //
  // Unless a particular reaction channel is represented by a data file given
  // in this list, it will not be included in the MARLEY simulation.
  reactions: [ "../react/ve40ArCC_Bhattacharya2009.react" ],

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
  //   Fermi-Dirac spectrum                   "fd", "fermi-dirac",
  //                                          "fermi_dirac"
  //
  //
  //   Pinched Fermi-Dirac spectrum           "bf", "beta", "beta-fit"
  //   with pinching parameter beta
  //   (see, e.g.,
  //   http://arxiv.org/abs/1511.00806)
  //
  //   Monoenergetic source                   "mono", "monoenergetic"
  //
  //   Muon decay at rest source              "dar", "decay-at-rest"
  //   (ve and vu only)
  //
  //   User-defined neutrino spectrum         "hist", "histogram"
  //   histogram
  //
  //   User-defined neutrino spectrum         "grid"
  //   probability density function
  //   described by an interpolating
  //   function on a set of grid points
  //
  //   ROOT neutrino spectrum histogram       "th1"
  //   (TH1)
  //
  //   ROOT neutrino spectrum probability     "tgraph"
  //   density function (TGraph)
  //
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
  //    eta: 0             // Pinching parameter (dimensionless, default 0)
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
  //  MUON DECAY AT REST
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
  // In this example configuration file, we've chosen a Fermi-Dirac source.
  //
  source: {
   type: "fermi-dirac",
   neutrino: "ve",       // The source produces electron neutrinos
   Emin: 0,              // Minimum neutrino energy (MeV)
   Emax: 60,             // Maximum neutrino energy (MeV)
   temperature: 3.5,     // Temperature (MeV)
   eta: 0                // Pinching parameter (dimensionless, default 0)
  },

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
  // In this example, incident neutrinos travel in the +z direction.
  direction: { x: 0.0,
               y: 0.0,
               z: 1.0
             },

  // LOGGER CONFIGURATION (optional)
  //
  // The "log" JSON array contains a list of JSON objects representing
  // one or more output streams that MARLEY should use to output diagnostic
  // messages while generating events.
  //
  // Each entry is a JSON object with the following keys:
  //
  //   - file: The name of a file to receive logger output. If the value
  //           is set to "stdout" or "stderr", then the corresponding stream
  //           will be used instead of a file on disk.
  //
  //   - level: The logging level that should be used as a threshold for
  //            writing to the stream. Valid values (in order of increasing
  //            severity) are "debug", "info", "warning", "error", and
  //            "disabled". Only  messages that are at least as severe as the
  //            given logging level will be written to the stream. If this key
  //            is omitted, a value of "info" will be assumed. A value of
  //            "disabled" suppresses all output from MARLEY to the stream.
  //
  //   - overwrite: Boolean value indicating whether any previously existing
  //                content in this stream should be erased by MARLEY (true)
  //                or not (false). This key is ignored when using stdout or
  //                stderr. If this key is omitted, a value of false is
  //                assumed.
  //
  // If the "log" key is omitted, MARLEY assumes a value of
  // [ { file: "stdout", level: "info" } ]
  //
  log: [ { file: "stdout", level: "info" },
         { file: "marley.log", level: "info", overwrite: true } ],

  // EXECUTABLE SETTINGS (optional)
  //
  // The entries within the executable_settings JSON object are used to
  // control the marley command-line executable. They are ignored if
  // the configuration file is used to initialize MARLEY outside of that
  // context (e.g., within a Geant4 application that links to the MARLEY
  // shared libraries).
  executable_settings: {

    // EVENT COUNT (optional)
    //
    // Generate 1e4 events (the JSON parser expects this entry to be an
    // integer literal, so scientific notation is not currently allowed for the
    // value).
    //
    // If this key is omitted, a value of 1000 will be assumed.
    events: 10000,


    // EVENT OUTPUT (optional)
    //
    // The "output" JSON array contains a list of JSON objects representing
    // one or more output streams that will receive the events generated
    // by the marley command-line executable
    //
    // Each entry is a JSON object with the following keys:
    //
    //   - file: The name of a file that will store the generated events.
    //           Streaming events to stdout or stderr is not currently
    //           supported.
    //
    //   - format: The format to use when storing the events in the file.
    //             Valid values are "ascii", "hepevt", "json", and "root".
    //             Details about the format options are given below.
    //
    //   - mode: The file I/O mode to use when writing to this file. For
    //           the "ascii" and "hepevt" formats, valid values are
    //           "overwrite" (erase any previously existing file contents)
    //           and "append" (continue output immediately after any
    //           existing file contents). For the "json" and "root" formats,
    //           valid values are "overwrite" and "resume". If the "resume"
    //           mode is chosen, the generator will restore its previous state
    //           from an incomplete run (e.g., a run that was interrupted
    //           by the user via ctrl+C) that was saved to the output file
    //           and continue from where it left off.
    //
    //   - force: Boolean value used only for the "overwrite" mode. If
    //            it is true, the marley executable will not prompt the
    //            user before overwriting existing data. If this key
    //            is omitted, a value of false is assumed.
    //
    //   - indent: Integer value used only for the "json" format. Gives the
    //             number of spaces that should be used as a tab stop when
    //             pretty-printing the JSON output. A value of 0 will result
    //             in the most compact output file, while nonzero values
    //             will be more human-readable. If this key is omitted, a
    //             value of 1 is assumed.
    //
    // The allowed output file formats are
    //
    //   - "ascii": The native format for MARLEY events. Files written in
    //              this format may be read and written by the marley::Event
    //              class using the >> and << stream operators, respectively.
    //
    //   - "hepevt": Each event is described by a HEPEVT format record.
    //               See section 3.1 of the StdHep manual
    //               (http://tinyurl.com/StdHepManual) for a description
    //               of the HEPEVT format.
    //
    //   - "json": The events are stored as an array of JSON objects. The
    //             state of the generator is also written to the file when
    //             execution terminates (either because the desired number
    //             of events has been reached or because the user has
    //             interrupted the run via ctrl+C). The function
    //             marley::Event::to_json() controls the output format.
    //
    //   - "root": Stores the generated marley::Event objects in a ROOT
    //             TTree. This format is only available if MARLEY has been
    //             built with ROOT support.
    //
    // If this key is omitted, then a value of
    // [ { file: "events.json", format: "json", mode: "overwrite" } ]
    // is assumed.
    //
    output: [ { file: "events.ascii", format: "ascii", mode: "overwrite" } ],
  },

} // A closing curly brace should appear at the end of the file
