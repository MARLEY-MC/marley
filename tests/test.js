{
  //energy_pdf_max: 0.000001,

  // Reaction matrix element files
  reactions: [ "ve40ArCC_Liu1998.react" ],

  // Neutrino source specification
  source: {
    type: "dar",
    neutrino: "ve",
  },

  //// Neutrino source specification
  //source: {
  //  type: "fd",
  //  neutrino: "ve",
  //  Emin: 0., 
  //  Emax: 60.,
  //  temperature: 4.5,
  //  eta: 10.,
  //},

  // Incident neutrino direction 3-vector
  direction: { x: 0.0,
               y: 0.0,
               z: 1.0
             },

  // Logging configuration
  log: [ { file: "stdout", level: "info" }, ],

  // Settings for marley command-line executable
  executable_settings: {

    events: 100000, // The number of events to generate

    // Event output configuration
    output: [
      { file: "test.root", format: "root", mode: "overwrite" },
      { file: "test.ascii", format: "ascii", mode: "overwrite" },
      { file: "test.hepevt", format: "hepevt", mode: "overwrite" },
      { file: "test.json", format: "json", mode: "overwrite" }
    ],
  },
}
