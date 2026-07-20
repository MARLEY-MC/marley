// Use this example job configuration file as a starting point for your own
// files.
{
  seed: 123456, // Random number seed (omit to use time since Unix epoch)

  // Pure 40Ar target
  target: {
    nuclides: [ 1000180400 ],
    atom_fractions: [ 1.0 ],
  },

  // Simulate ve scattering on 40Ar:
  reactions: [
    // Charged-current interactions on the nucleus using the
    // recommended MARLEY v2 cross-section models for both
    // continuum (HF-CRPA) and discrete (Bhattacharya2009)
    // nuclear transitions
    "ve40ArCC_HF-CRPA.react",
    "ve40ArCC_Bhattacharya2009-Discrete.react",
    // Elastic scattering on atomic electrons
    "ES.react"
  ],

  // Neutrino source specification
  source: {
    type: "fermi-dirac",
    neutrino: "ve",      // The source produces electron neutrinos
    Emin: 0,             // Minimum neutrino energy (MeV)
    Emax: 60,            // Maximum neutrino energy (MeV)
    temperature: 3.5,    // Temperature (MeV)
    eta: 4               // Pinching parameter (dimensionless, default 0)
  },

  // Incident neutrino direction 3-vector
  direction: { x: 0.0, y: 0.0, z: 1.0 },

  // Settings for the "marley generate" command
  generate: {

    // The number of events to generate
    events: 10000,

    // Event output configuration
    output: [ { file: "events.hepmc3", format: "ascii", mode: "overwrite" } ],
  },
}
