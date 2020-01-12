# Example ROOT macros

This directory contains several example ROOT macros that can be used to plot
truth information from MARLEY events.  To use these macros, you will need to
build MARLEY against an installation of ROOT. This will be handled by the
standard Makefile (${MARLEY}/build/Makefile) automatically if the root-config
executable is visible on the system PATH.

Once MARLEY has been linked to ROOT successfully, running the mroot script
(located in ${MARLEY}/build/ following a successful build) will start the ROOT
interpreter and load class dictionaries needed to interpret marley::Event
objects. All of the macros below are designed to be run from within an mroot
session. The input file may contain MARLEY events in any of the four standard
output formats (ASCII, HEPEVT, JSON, and ROOT).

Brief descriptions of each of the example macros, together with the syntax
needed to run them within the mroot interpreter, are given below.

cos_plot.C
-----------

Plots the flux-averaged differential cross section with respect to the
scattering cosine of the outgoing electron.

Usage:

  .x cos_plot.C("/path/to/event/file")

Ex.C
----

Plots a histogram of excitation energies accessed by the initial two-two
scatter on the target nucleus.

Usage:

  .x Ex.C("/path/to/event/file")

fp_spect.C
----------

Plots a histogram of kinetic energies for the final particle type with
a given Particle Data Group (PDG) code.

Usage:

  .x fp_spect.C("/path/to/event/file", PDG code)


nu_spect.C
----------

Plots the flux-averaged differential cross section with respect to the true
neutrino energy.

Usage:

  .x nu_spect.C("/path/to/event/file")

print_events.C
--------------

Prints a human-readable description of each event in the input file.

Usage:

  .x print_events.C("/path/to/event/file")

reco_spectrum.C
---------------

Plots the flux-averaged differential cross section with respect to the true
neutrino energy (black) and two reconstructed neutrino energies. The "reco 1"
neutrino energy (blue) is calculated from the sum of the true kinetic energies
of all final-state particles excluding neutrons. The "reco 2" neutrino energy
(red) includes only the kinetic energy of the final-state electron. Both
"reco 1" and "reco 2" also include the energy needed for a transition from the
ground state of the target 40Ar nucleus to the ground state of the daughter 40K
nucleus.

Usage:

  .x reco_spectrum.C("/path/to/event/file")
