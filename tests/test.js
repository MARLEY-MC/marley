{
  // Reaction matrix element files
  reactions: [ "CEvNS40Ar.react" ],

  // Neutrino source specification
  source: {
    type: "dar",
    neutrino: "vubar",
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
}
