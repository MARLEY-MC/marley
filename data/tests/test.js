{
  // Reaction matrix element files
  reactions: [ "CEvNS40Ar.react" ],

  // Use the allowed approximation for CEvNS form factors.
  // The analytic truth in src/tests/sampling.cc assumes this simplified model
  // (gV2=1, no momentum-transfer-dependent form factors), so the generator
  // configuration must match to make the chi-squared test meaningful.
  form_factors: "AA",

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
