{
  // Reaction matrix element files
  reactions: [ "CEvNS40Ar.react" ],

  // Neutrino source specification
  source: {
    type: "dar",
    neutrino: "vubar",
  },

  // Incident neutrino direction 3-vector
  direction: { x: 0.0,
               y: 0.0,
               z: 1.0
             },

  // Logging configuration
  log: [ { file: "stdout", level: "info" }, ],
}
