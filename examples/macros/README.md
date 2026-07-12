# Example ROOT macros

This directory contains several example ROOT macros that can be used to plot
truth information from MARLEY events.

These macros consume the ROOT summary file produced by the `marley summarize`
command. To use them, first run

```
marley summarize -o <output.root> <input_events_1> ...
```

where the input event file(s) may be in any of the HepMC3 formats supported by
MARLEY (ASCII or ROOT). Once the summary ROOT file has been created, open it in
ROOT and run any of the macros below:

```
root
root [0] .x cos_plot.C("output.root")
```

Brief descriptions of each macro and its syntax are given below.

cos_plot.C
-----------

Plots the flux-averaged differential cross section with respect to the
scattering cosine of the outgoing lepton.

Usage:

  .x cos_plot.C("/path/to/summary.root")

Ex.C
----

Plots a histogram of excitation energies accessed by the initial two-to-two
reaction on the target nucleus.

Usage:

  .x Ex.C("/path/to/summary.root")

fp_spect.C
----------

Plots a histogram of kinetic energies for the final particle type with
a given Particle Data Group (PDG) code.

Usage:

  .x fp_spect.C("/path/to/summary.root", PDG code)


nu_spect.C
----------

Plots the flux-averaged differential cross section with respect to the true
neutrino energy.

Usage:

  .x nu_spect.C("/path/to/summary.root")

reco_spectrum.C
---------------

This macro is intended for use with a sample of events involving CC electron
neutrino scattering on 40Ar only.

It plots the flux-averaged differential cross section with respect to the true
neutrino energy (black) and two reconstructed neutrino energies. The "reco 1"
neutrino energy (blue) is calculated from the sum of the true kinetic energies
of all final-state particles excluding neutrons. The "reco 2" neutrino energy
(red) includes only the kinetic energy of the final-state electron. Both
"reco 1" and "reco 2" also include the energy needed for a transition from the
ground state of the target 40Ar nucleus to the ground state of the daughter 40K
nucleus.

Usage:

  .x reco_spectrum.C("/path/to/summary.root")
