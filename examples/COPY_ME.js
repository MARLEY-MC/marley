// Use this example configuration file as a starting point for your own files.
{
  seed: 123456, // Random number seed (omit to use time since Unix epoch)

  // Nuclear structure data files
  structure: [ "../structure/z019",
               "../structure/z018",
               "../structure/z017", ],

  // Reaction matrix element files
  reactions: [ "../react/ve40ArCC_Bhattacharya2009.react" ],

  // Neutrino source specification
  source: {
   type: "fermi-dirac",
   neutrino: "ve",       // The source produces electron neutrinos
   Emin: 0,              // Minimum neutrino energy (MeV)
   Emax: 60,             // Maximum neutrino energy (MeV)
   temperature: 3.5,     // Temperature (MeV)
   eta: 0                // Pinching parameter (dimensionless, default 0)
  },

  // Incident neutrino direction 3-vector
  direction: { x: 0.0,
               y: 0.0,
               z: 1.0
             },

  // Logging configuration
  log: [ { file: "stdout", level: "info" },
         { file: "marley.log", level: "info", overwrite: true } ],

  // Settings for marley command-line executable
  executable_settings: {

    events: 10000, // The number of events to generate

    // Event output configuration
    output: [ { file: "events.ascii", format: "ascii", mode: "overwrite" } ],
  },
}
