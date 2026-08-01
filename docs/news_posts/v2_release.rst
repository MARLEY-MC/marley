MARLEY v2
=========

.. feed-entry::
   :date: 2026-08-01 15:11
   :author: Steven Gardiner

I am pleased to announce `version 2.0.0
<https://github.com/MARLEY-MC/marley/releases/tag/v2.0.0>`__ of MARLEY, the
first official update to the code since November 2023. This major release
implements an improved model of the :math:`\nu_e\text{-}^{40}\mathrm{Ar}`
charged-current cross section that is described in |arXiv|. The previous MARLEY
output format has been entirely replaced by the `HepMC3
<https://arxiv.org/abs/1912.08005>`__ event record. The representation of
neutrino interactions within HepMC3 data structures follows the recommendations
of version 1.0.0 of the `NuHepMC
<https://doi.org/10.21468/SciPostPhysCodeb.57>`__ standard.

.. cut::

Physics updates
---------------

The headline change in v2.0.0 is an overhaul of MARLEY's neutrino cross-section
calculation that removes the rough *allowed approximation* (AA) treatment used
in the v1 release series. In MARLEY v2, neutrino-induced transitions to
high-lying nuclear excitation energies are modeled using a Hartree-Fock
Continuum Random Phase Approximation (HF-CRPA) approach, including forbidden
nuclear matrix elements. The data-driven strategy from v1 is retained for
transitions to discrete nuclear levels, but approximate corrections for
momentum-transfer dependence are introduced to avoid taking the :math:`q \to 0`
limit from the AA that was used previously. The updated hybrid model of the
charged-current interaction is discussed at length in |arXiv|.

Hauser-Feshbach statistical de-excitation calculations in MARLEY v2 still use
the 46-parameter nuclear optical potential developed for the 2003
Koning-Delaroche `global fit
<https://doi.org/10.1016/S0375-9474(02)01321-0>`__. However, the default
parameter values have now been updated to match the "federal" `KDUQ
<https://doi.org/10.1103/PhysRevC.107.014602>`__ re-evaluation published in
2023. All parameters used in the nuclear optical model are now configurable via
the JSON job configuration file used to run MARLEY.

Output events now contain the full nuclear de-excitation history rather than
merely vectors of the initial- and final-state particles. Each binary decay
vertex in the de-excitation cascade is explicitly recorded. Metadata fields in
the event record are used to store properties of the intermediate nuclear
states (excitation energy, spin, parity) as well as the total and partial decay
widths associated with each sampled transition.

Rather than treating nuclear de-excitations as instantaneous (as was done in
v1), MARLEY v2 adds tracking of particle emission times to the simulation of
continuum decays and (where half-life data are available) γ-ray transitions
between discrete levels. This functionality had to be introduced externally to
support an exciting recent analysis by DEAP-3600 that reported first evidence
for charged-current absorption of :math:`^{8}\mathrm{B}` solar neutrinos on an
argon target: |deap|. It is now available "out of the box" for all users.

New code infrastructure for event reweighting has been added that enables
evaluation of systematic uncertainties without re-running the simulation. A
full uncertainty budget for MARLEY predictions has not yet been developed, but
weight calculators for two classes of uncertainties (related to discrete
transition matrix elements and the KDUQ nuclear optical potential) are
available for the first time in version 2.0.0.

Technical updates
-----------------

The MARLEY-specific output format used in v1 has been deprecated in favor of
`HepMC3 <https://arxiv.org/abs/1912.08005>`__, a highly-flexible event record
definition that is widely used for collider physics simulations. Recent
community discussions about neutrino event generator data formats led to
creation of `NuHepMC <https://doi.org/10.21468/SciPostPhysCodeb.57>`__, a set
of conventions for describing simulated neutrino reactions within the HepMC3
framework. The updated MARLEY event format follows version 1.0.0 of the NuHepMC
standard. **This is a breaking change** from v1 that will typically require
updates to later stages of experimental simulation chains and downstream
analysis code. The `interpreting the output <../interpret_output.html>`__
webpage provides details about the new format and MARLEY's built-in
capabilities for conversions to other representations of the simulation
results.

The ``marley`` executable now provides a unified command-line interface that
accepts multiple commands: ``generate``, ``decay``, ``convert``, ``print``,
``reweight``, ``summarize``, ``xsec``, ``help``, and ``version``. These
commands replace the previous system of several different executables and
helper scripts. For backwards compatibility, ``marley CONFIG_FILE`` still works
as shorthand for ``marley generate CONFIG_FILE``. A new `command reference
<../commands.html>`__ webpage provides documentation for all of the new
commands.

MARLEY can now be built using `CMake <https://cmake.org>`__ in addition to
plain GNU Make. For convenience, the top-level ``Makefile`` will delegate to
CMake automatically if it is available and will fall back to a standalone GNU
Make recipe otherwise; users who prefer to invoke CMake directly are free to
follow the usual CMake build procedure instead. Note that a C++17-compliant
compiler (GCC ≥ 9.1.0 or Clang ≥ 9.0.0) is now required to build MARLEY.

The JSON-like job configuration file format has been extended to support an
``#include:"file"`` directive, allowing MARLEY's runtime configuration to be
split across multiple files. Error messages produced by the JSON parser now
report the file name and the line/column position at which a problem was
encountered, making configuration mistakes easier to diagnose and fix.

Outlook
-------

Thanks to a 2025 Early Career `award`_ from the U.S. Department of Energy,
MARLEY has recently returned to a state of active development. Further details
about physics model improvements in MARLEY v2 can be found in my
Neutrino 2026 `poster`_ presentation and `this talk`_ from the July 2026
APS DPF `meeting <https://indico.fnal.gov/event/72820/>`__.

For the v2.0.0 release, I thank all the co-authors of |arXiv| for their
scientific collaboration. I am also grateful for help with the implementation
of new MARLEY features from Pablo Barham Alzás and Luca Abu El-Haj
(cross-section improvements, optical potential reweighting) as well as Everardo
Granados Vázquez (γ-ray emission timing).

The MARLEY code remains enthusiastically open-source. `Feedback
<../contact.html>`__ and `contributions <../dev_docs.html>`__ from the
community are both valued and welcomed. I look forward with excitement as you
give MARLEY v2.0.0 a try!

.. _`arXiv:2604.26801 [hep-ph]`: https://arxiv.org/abs/2604.26801

.. |arXiv| replace:: `arXiv:2604.26801 [hep-ph]`_

.. _`arXiv:2605.12769 [hep-ex]`: https://arxiv.org/abs/2605.12769

.. |deap| replace:: `arXiv:2605.12769 [hep-ex]`_

.. _`award`: https://news.fnal.gov/2026/04/steven-gardiner-receives-early-career-award-to-advance-low-energy-neutrino-research-at-dune/

.. _`poster`: https://indico.global/event/15740/contributions/147631/attachments/72141/139992/Neutrino-2026-Poster-final.pdf

.. _`this talk`: https://indico.fnal.gov/event/72820/contributions/341239/attachments/199988/278787/marley-v2-aps-dpf.pdf
