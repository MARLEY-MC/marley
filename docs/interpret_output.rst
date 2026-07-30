=======================
Interpreting the output
=======================

This page provides detailed guidance for parsing and analyzing output files
produced by the ``marley`` command-line executable. It begins with a
description of the standard *PDG codes* used to identify particle types in all
MARLEY simulation results. This is followed by documentation for each available
output file format.

.. _pdg-codes:

PDG codes
^^^^^^^^^

The `Particle Data Group <https://pdg.lbl.gov/index.html>`__ (PDG) has defined
a standard numbering scheme for representing particle species in Monte Carlo
event generators. Each kind of particle is assigned a unique positive integer
as an identifier. The corresponding antiparticle is assigned a negative integer
with the same absolute value. A full description of the numbering scheme is
available `here <https://pdg.lbl.gov/current/mc-particle-id>`__.

Like nearly all modern event generators used in particle physics, MARLEY adopts
the integer *PDG codes* for particle identification and uses them both
internally and in output files. For convenience, a table of the PDG codes most
relevant for MARLEY is given below.

.. list-table::
   :header-rows: 1
   :class: pdg-table

   * - PDG code
     - Particle
   * - 11
     - :math:`e^-`
   * - 12
     - :math:`\nu_e`
   * - 13
     - :math:`\mu^-`
   * - 14
     - :math:`\nu_\mu`
   * - 15
     - :math:`\tau^-`
   * - 16
     - :math:`\nu_\tau`
   * - 22
     - :math:`\gamma`
   * - 2112
     - :math:`n`
   * - 2212
     - :math:`p`
   * - 1000010020
     - :math:`d`
   * - 1000010030
     - :math:`t`
   * - 1000020030
     - :math:`h`
   * - 1000020040
     - :math:`\alpha`

In general, a PDG code of the form 100ZZZAAA0 represents a nuclide with proton
number Z and mass number A. For example, :math:`^{40}\mathrm{Ar}` is
represented by the PDG code 1000180400.

The HepMC3 event graph
^^^^^^^^^^^^^^^^^^^^^^

The data structures defined by the `HepMC3
<https://arxiv.org/abs/1912.08005>`__ event record `library
<https://gitlab.cern.ch/hepmc/HepMC3>`__ are used to implement the canonical
representation of physics events in MARLEY v2.0.0 and later. To facilitate
interoperability with other software tools used in neutrino physics, MARLEY
also implements version 1.0.0 of the `NuHepMC
<https://doi.org/10.21468/SciPostPhysCodeb.57>`__ standard, which defines a set
of conventions for representing neutrino interaction events within the HepMC3
infrastructure in a generator-agnostic way.

A key concept for reading MARLEY output files is the HepMC3 *event graph*. Each
event is represented as a directed graph in which **particles** are the edges
and **vertices** are the nodes. A vertex groups one or more incoming particles
with one or more outgoing particles and typically represents a physical
interaction. Every particle in a HepMC3 event is assigned a unique positive
integer ID, and every vertex is assigned a unique negative integer ID. A
particle that has no production vertex (i.e., an initial-state particle) is
indicated by a production vertex ID of zero. Note that these IDs are used to
index individual particles in the event graph; they are distinct from the
:ref:`PDG codes <pdg-codes>` mentioned above that are used to differentiate
particle species.

In a typical MARLEY neutrino-nucleus scattering event, the HepMC3 graph
contains two kinds of vertices. The primary vertex represents the hard 2-to-2
interaction that produces the final-state lepton and outgoing nucleus. One or
more de-excitation vertices may also be present in which the outgoing nucleus
emits γ-rays, nucleons, or light complex fragments. Vertices representing
particle emissions from the unbound continuum and from discrete nuclear levels
are assigned distinct status codes, and a separate vertex is added to the event
record for each binary decay step in the de-excitation cascade until the
nuclear ground state is reached.

Output file formats
^^^^^^^^^^^^^^^^^^^

The ``generate``, ``reweight``, and ``decay`` commands accepted by the
``marley`` executable produce output files containing physics events stored
within HepMC3 data structures. The generic description of a HepMC3 event graph
can be serialized in multiple ways, and MARLEY currently supports both an ASCII
text representation and (when built with ROOT enabled) a binary ROOT-based
representation as official, full-featured output formats. The ``convert``
command can translate between these equivalent formats and also provide one-way
conversions to some of the deprecated event formats used in MARLEY v1.2.1 and
earlier. The ``summarize`` command creates a simplified "flat" ROOT `TTree
<https://root.cern/doc/master/group__tree.html>`__ that stores a subset of the
full event information for more convenient analysis. It replaces the
functionality provided by the ``marsum`` executable that existed in the MARLEY
v1 release series. Descriptions of the various output formats produced by these
commands are given below.


ASCII HepMC3
------------

A file in the ASCII HepMC3 format begins with a two-line header::

  HepMC::Version 3.02.07
  HepMC::Asciiv3-START_EVENT_LISTING

followed by a single run-info block and then one event block per simulated
event. The file ends with the line::

  HepMC::Asciiv3-END_EVENT_LISTING

Run-info block
~~~~~~~~~~~~~~

The run-info block appears once at the start of the file and records information
about the MARLEY run as a whole. Each line begins with a single-character tag.

``W`` — weight names
  A backslash-pipe (``\|``) separated list of the weight names defined for
  this run. MARLEY always declares at least ``CV`` (central value). If weight
  calculators were configured and used during a ``generate`` or ``reweight``
  job, then additional weight names may appear later on this line::

    W CV

``T`` — tool identification
  The generator name (MARLEY), version, and git commit hash, separated by
  backslash-pipes.

``A`` — run-level attributes
  Each ``A`` line stores one named attribute as ``A <name> <value>``.
  MARLEY writes the following run-level attributes:

  ``MARLEY.JSONconfig``
    The complete job configuration used to generate the events, serialized
    as a JSON string. This makes the file self-describing: all settings
    (reaction files, neutrino source, random seed, output paths, etc.) are
    recorded.

  ``MARLEY.RNGseed``
    The integer seed that was used to initialize the random number generator.

  ``MARLEY.ReweightConfig.<n>`` (reweighted files only)
    The reweighting configuration applied by ``marley reweight``, serialized as
    JSON. One attribute (indexed by the integer ``<n>``) is written per
    reweighting pass applied.

  ``NuHepMC.Version.{Major,Minor,Patch}``
    The version of the NuHepMC standard implemented (currently 1, 0, 0).

  ``NuHepMC.Conventions``
    A space-separated list of the NuHepMC convention labels that MARLEY
    follows. For v2.0.0 output this is always::

      G.C.2 G.C.3 E.C.1 E.C.2 E.C.3

    These labels have the following meaning:

    - **G.C.2**: The flux-averaged total cross section :math:`\langle \sigma
      \rangle` is known before the run begins and is stored once in the run-info
      block (see ``NuHepMC.FluxAveragedTotalCrossSection`` below), rather than
      being updated in each event as the generator runs.
    - **G.C.3**: Citation metadata for MARLEY's physics models is embedded in
      the run-info block (see ``NuHepMC.Citations.*`` below).
    - **E.C.1**: Process IDs follow the NuHepMC recommended identifier ranges
      (100–199 for low-energy nuclear scattering, etc.; see the process ID
      table in the ``marley summarize`` section below).
    - **E.C.2**: Each event records the total interaction cross section for
      the projectile at its sampled energy (``tot_xs``).
    - **E.C.3**: Each event records the partial cross section for the selected
      primary interaction process (``proc_xs``).

  ``NuHepMC.FluxAveragedTotalCrossSection``
    The flux-averaged total cross section :math:`\langle \sigma \rangle` for the
    run, in picobarns per target atom. This is the quantity needed to convert a
    distribution of simulated events into a cross-section prediction.

  ``NuHepMC.Citations.Generator.{DOI,InspireHEP,arXiv}``
    Space-separated lists of the DOIs, InspireHEP keys, and arXiv identifiers
    for the papers describing MARLEY's physics models. These should be cited
    whenever MARLEY output is used in a publication.

  ``NuHepMC.ProcessIDs`` and ``NuHepMC.ProcessInfo[<ID>].{Name,Description}``
    The complete list of process IDs that may appear in events from this run,
    together with a name and description for each. Reading these entries
    directly from the file is the definitive way to interpret the
    ``signal_process_id`` event attribute.

  ``NuHepMC.VertexStatusIDs`` and ``NuHepMC.VertexStatusInfo[<ID>].{Name,Description}``
    The vertex status codes used in the file and their meanings. MARLEY
    uses three codes: 1 (``Primary``, the primary interaction vertex), 22
    (``HFDecay``, a continuum de-excitation step), and 23
    (``GammaDecay``, a de-excitation step from a bound nuclear energy level).

  ``NuHepMC.ParticleStatusIDs`` and ``NuHepMC.ParticleStatusInfo[<ID>].{Name,Description}``
    The particle status codes used in the file and their meanings. MARLEY
    uses five codes: 1 (``Final-state``), 4 (``Projectile``), 20
    (``Target``), 27 (``UndecayedRemnant``, a nuclear residue before
    de-excitation begins), and 28 (``IntermediateRemnant``, a nucleus
    undergoing a de-excitation cascade).

  ``NuHepMC.Units.CrossSection.Unit`` and ``NuHepMC.Units.CrossSection.TargetScale``
    The units for all cross section values in the file. MARLEY always uses
    picobarns (``pb``) per target atom (``PerAtom``).

  ``NuHepMC.AdditionalParticleNumbers`` (and ``NuHepMC.AdditionalParticleNumbers[<PDG>].{Name,Description}``)
    Particle codes used in the file that are not in the standard PDG
    numbering scheme. MARLEY declares PDG code 0, the ``Absent`` dummy
    projectile used in ``marley decay`` output (where there is no real
    incoming beam particle).

Per-event blocks
~~~~~~~~~~~~~~~~

The first event block begins immediately after the run-info block on a line
starting with an ``E`` character. An event block ends at the next ``E`` line or
the end-of-file marker. The lines within an event block have the following
structure.

``E`` — event header::

  E <event_number> <vertex_count> <particle_count>

For example, ``E 1 5 12`` opens event number 1, which contains 5 vertices
and 12 particles.

``U`` — units::

  U MEV CM

MARLEY always uses MeV for energies and cm for positions.

``W`` — event weights::

  W 1.0000000000000000000000e+00

One space-delimited weight value for each name declared in the run-info ``W``
line. For unweighted output from ``marley generate``, only the ``CV`` (for
"central-value") weight appears and is always exactly one. In cases where
additional event weight have been computed, their numerical values appear in the
same order as the weight names listed in the run information.

``A`` — per-event and per-particle attributes
  Each attribute line has the form ``A <scope> <name> <value>``, where the
  scope identifies the object to which the attribute belongs:

  - **0**: an event-level attribute.
  - **Positive integer**: an attribute attached to the particle with this ID.
  - **Negative integer**: an attribute attached to the vertex with this ID.

  *Event-level attributes (scope 0):*

  ``lab_pos``
    The interaction position in the lab frame as three space-separated values
    (x, y, and z; all in cm). Zeros are always given unless MARLEY has been
    interfaced with a separate detetor simulation.

  ``signal_process_id``
    An integer code identifying the primary interaction type for this event. The
    meaning of each signal_process_id value is defined in the run-info block.

  ``tot_xs``
    The total interaction cross section summed over all active processes,
    evaluated at the projectile's sampled energy for the current event (pb per
    atom).

  ``proc_xs``
    The partial cross section for the specific process selected as the primary
    interaction. This is evaluated at the projectile's sampled energy for the
    current event (pb per atom).

  ``MARLEY.GeneratorState``
    A serialized snapshot of the random number generator state. This attribute
    is written only by the ``generate`` command, and it appears only on the last
    event in the file. It is saved upon normal job completion or early
    termination (unhandled exception or user interrupt via ctrl+C). Its purpose
    is to support the ``resume`` output mode, which allows an interrupted
    ``generate`` job to continue exactly where it left off. For information
    about how to configure the ``resume`` output mode, see the example
    ``generate`` job configuration file (``examples/config/annotated.js``).

  *Particle attributes (positive scope = particle ID):*

  ``Ex``
    Nuclear excitation energy (MeV). This attribute is stored for particles
    representing an atomic nucleus. A value of 0 indicates the ground state.

  ``twoJ``
    Two times the nuclear spin quantum number J. Doubling the value of J allows
    half-integer spins to be represented as integers. This attribute is stored
    for particles representing an atomic nucleus.

  ``parity``
    Intrinsic parity (±1). This attribute is stored for particles representing
    an atomic nucleus.

  ``charge``
    Net electric charge of the particle (integer, in units of the elementary
    charge). Written on atomic/ionic particles: the target atom and the outgoing
    nucleus at each stage of the de-excitation cascade. This records the charge
    of the atom or ion as a whole rather than of the bare nucleus, which is
    already encoded in the PDG code. When the particle's charge can be
    unambiguously determined by PDG code (as it can for elementary particles),
    this attribute is not stored.

``P`` — particle lines
  Each particle in the event is described by a ``P`` line::

    P <id> <prod_vtx_id> <pdg> <px> <py> <pz> <E> <mass> <status>

  ``id`` is the particle's unique positive integer identifier within the event.
  ``prod_vtx_id`` is the ID of the vertex that produced this particle (negative
  integer), or 0 for initial-state particles that have no production vertex.
  ``pdg`` is the particle's PDG code, ``px``, ``py``, ``pz``, ``E``, and
  ``mass`` are the 3-momentum components, total energy, and mass, all in MeV.
  ``status`` is one of the particle status codes defined in the run-info block.

``V`` — vertex lines
  Each vertex is described by a ``V`` line::

    V <id> <status> [<parent_particle_ids>]

  or, for vertices with a non-zero spacetime position::

    V <id> <status> [<parent_particle_ids>] @ <x> <y> <z> <t>

  ``id`` is the unique negative integer identifier for the vertex. ``status`` is
  one of the vertex status codes defined in the run-info block. The bracketed
  list gives the IDs of the particles entering this vertex (i.e., particles for
  which this is their end vertex). The optional ``@ x y z t`` suffix gives the
  spacetime position in cm and cm/c; it is present for de-excitation vertices
  that have a non-zero time delay and omitted otherwise.

Worked example
~~~~~~~~~~~~~~

The following ASCII-format HepMC3 event represents a charged-current primary
interaction :math:`\nu_e + \, ^{40}\mathrm{Ar} \to e^{-} + \, ^{40}\mathrm{K}^*`
followed by a chain of four γ-ray emissions from nuclear de-excitations. Momenta
and energies in the ``P`` lines are rounded to three decimal places for
readability; an actual MARLEY output file would use enough digits to preserve
full double-precision for all floating-point numbers.

::

  E 1 5 12
  U MEV CM
  W 1.0000000000000000000000e+00
  A 4 Ex 4.3837
  A 6 Ex 2.28987
  A 8 Ex 1.64364
  A 10 Ex 0.0298299
  A 12 Ex 0
  A -5 GammaBranchingRatio 1
  A -4 GammaBranchingRatio 0.803859
  A -3 GammaBranchingRatio 0.562746
  A -2 GammaBranchingRatio 0.757576
  A -5 TotalWidth 1.07350060716826e-13
  A -4 TotalWidth 1.35785047037652e-15
  A -3 TotalWidth 5.4968404583917e-09
  A 2 charge 0
  A 4 charge 1
  A 6 charge 1
  A 8 charge 1
  A 10 charge 1
  A 12 charge 1
  A 0 lab_pos 0.000000 0.000000 0.000000
  A 4 parity 1
  A 6 parity 1
  A 8 parity 1
  A 10 parity -1
  A 12 parity -1
  A 0 proc_xs 6.05325088059723e-05
  A 0 signal_process_id 100
  A 0 tot_xs 7.45272521695649e-05
  A 4 twoJ 0
  A 6 twoJ 2
  A 8 twoJ 0
  A 10 twoJ 6
  A 12 twoJ 8
  P 1 0 12 0.000 0.000 21.200 21.200 0.000 4
  P 2 0 1000180400 0.000 0.000 0.000 37224.7 37224.7 20
  V -1 1 [1,2]
  P 3 -1 11 -7.848 -7.396 11.564 15.821 0.511 1
  P 4 -1 1000190400 7.848 7.396 9.636 37230.1 37230.1 27
  V -2 23 [4]
  P 5 -2 22 0.134 1.298 -1.637 2.094 0.000 1
  P 6 -2 1000190400 7.713 6.098 11.273 37228.0 37228.0 28
  V -3 23 [6] @ 0.000 0.000 0.000 1.001e-03
  P 7 -3 22 0.373 -0.320 0.419 0.646 0.000 1
  P 8 -3 1000190400 7.340 6.419 10.854 37227.4 37227.4 28
  V -4 23 [8] @ 0.000 0.000 0.000 4.150e+03
  P 9 -4 22 -0.834 -0.559 1.264 1.614 0.000 1
  P 10 -4 1000190400 8.174 6.978 9.590 37225.7 37225.7 28
  V -5 23 [10] @ 0.000 0.000 0.000 4.294e+03
  P 11 -5 22 0.019 -0.021 0.010 0.030 0.000 1
  P 12 -5 1000190400 8.154 6.999 9.580 37225.7 37225.7 1

Walking through this event:

- **Particles 1 and 2** are the initial-state particles: an electron neutrino
  (PDG 12, status 4 = Projectile) with kinetic energy ≈ 21.2 MeV travelling
  along the z-axis, and a stationary :sup:`40`\ Ar nucleus (PDG 1000180400,
  status 20 = Target) with rest mass ≈ 37224.7 MeV. Their production vertex
  ID is 0, indicating they are initial-state particles with no production
  vertex in the graph.

- **Vertex −1** (status 1 = Primary) represents the primary interaction,
  with particles 1 and 2 as inputs. It produces two outgoing particles.

- **Particle 3** is the final-state electron (PDG 11, status 1 = Final-state,
  mass ≈ 0.511 MeV). **Particle 4** is the :math:`^{40}\mathrm{K}^*` nuclear
  residue (PDG 1000190400, status 27 = UndecayedRemnant) in an excited state
  with excitation energy :math:`E_x = 4.3837 \; \mathrm{MeV}`,
  spin :math:`J = 0` (``twoJ`` = 0), and positive parity. This discrete nuclear
  level is the isobaric analog of the :math:`^{40}\mathrm{Ar}` ground state.

- **Vertices −2 through −5** (all status 23 = GammaDecay) are successive
  de-excitation steps simulated using tabulated γ-ray branching ratios. Each
  takes the intermediate nuclear remnant (status 28) as input and produces one
  γ-ray (PDG 22, status 1) and a new remnant at a lower excitation energy.
  Vertices −3, −4, and −5 carry an ``@ … t`` suffix showing the time of each
  de-excitation step in cm/c; the first step (vertex −2) has no time suffix,
  meaning it is treated as instantaneous (due to an unknown nuclear level
  half-life). The ``Ex`` attributes on particles 4, 6, 8, 10, and 12 trace the
  excitation energy at each stage of the cascade: 4.3837 → 2.28987 → 1.64364 →
  0.0298299 → 0 MeV.

- **Particle 12** (status 1 = Final-state) is the :math:`^{40}\mathrm{K}`
  nucleus in its ground state after the cascade is complete.

The ``signal_process_id`` of 100 identifies this as a charged-current
interaction that populates a discrete nuclear energy level
(``vCC-discrete`` as defined in the run-info block).
The ``proc_xs`` and ``tot_xs`` attributes give the cross sections at this
event's neutrino energy: the process cross section ≈ 6.05 × 10\ :sup:`−5` pb
and the total cross section ≈ 7.45 × 10\ :sup:`−5` pb per target atom.

ROOT HepMC3
-----------

If MARLEY has been built with ROOT support (see the :doc:`getting_started` page
for build instructions), full HepMC3 events can also be stored in ROOT's
compressed binary format. The ROOT HepMC3 format is an equally full-featured
representation of the events, not merely a subset or summary.

MARLEY uses its own reader/writer classes (``OutputFileRoot`` and
``EventFileReader``) built around the ``HepMC3::GenEventData`` and
``HepMC3::GenRunInfoData`` plain-data structs from the official HepMC3 library.
It does not use HepMC3's ``WriterRoot`` or ``ReaderRoot`` classes.

A ROOT HepMC3 file produced by MARLEY contains:

- A ``TTree`` named ``MARLEY_event_tree`` with a single branch named ``event``
  holding one ``HepMC3::GenEventData`` object per event.

- A ``HepMC3::GenRunInfoData`` object named ``MARLEY_run_info``
  containing the run-level metadata (the same information as the run-info
  block in the ASCII HepMC3 format).

Opening a MARLEY ROOT HepMC3 file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After sourcing the ``setup_marley.sh`` environment script, one may open a
MARLEY ROOT file directly in a ROOT session::

  source setup_marley.sh
  root -l some_events.root

The ``libMARLEY`` library and its ROOT dictionary are loaded automatically by
ROOT's autoload mechanism. One can then browse the raw event data, for
example::

  MARLEY_event_tree->Scan("particles.pid:particles.status", "", "", 20)

The ``GenEventData`` struct fields available for browsing include
``particles`` (a vector of ``GenParticleData`` structs with fields
``pid``, ``momentum``, ``status``, etc.), ``vertices`` (a vector of
``GenVertexData`` structs), and ``event_number``, among others.

For most analysis purposes, the recommended approach is to use one of
MARLEY's built-in tools rather than navigating the HepMC3 data structures
directly:

- ``marley print`` — display events in a human-readable format.
- ``marley convert`` — translate between output formats.
- ``marley summarize`` — produce a flat ROOT ntuple (described below)
  for convenient histogram-level analysis.

Advanced users who need more complete access can write their own analysis code
that links against the HepMC3 library (``libHepMC3``) directly, using
``HepMC3::GenEvent::read_data()`` to reconstruct full event objects from the
stored ``GenEventData`` structs. This requires the HepMC3 headers and shared
library to be available at build time and runtime.

.. _summary-ntuple-format:

Summary ROOT TTree
------------------

The ``summarize`` command reads one or more MARLEY HepMC3 event files (in either
the ASCII or ROOT format) and writes a flat ROOT ntuple suitable for quick
analysis. This is a convenience projection of the full HepMC3 event record: it
removes many details of the event history while exposing the most
commonly-needed information. The summary ROOT TTree cannot be converted back
into the full HepMC3 event format.

Usage::

  marley summarize -o summary.root events.hepmc3

Multiple input files may be given; they are concatenated in the order specified
after checking for consistency of the configurations given in each run-info
block. The output file contains a single ``TTree`` named ``mst`` ("MARLEY
summary tree") with one row per event. A ``std::vector< std::string >`` object
called ``MARLEY_other_weight_names`` also appears in the file and stores the
names of any weights beyond the central-value weight (``CV``) that appeared in
the original HepMC3 events.

The example ROOT macros in ``examples/macros/`` all consume the ``mst`` tree
produced by ``marley summarize``; see ``examples/macros/README.md`` for details.
None of these example macros operate on raw ROOT HepMC3 files directly.

Branch listing
~~~~~~~~~~~~~~

.. |doubleType| raw:: html

   <i style="font-weight: normal;">(double)</i>

.. |intType| raw:: html

   <i style="font-weight: normal;">(int)</i>

.. |doubleVecType| raw:: html

   <i style="font-weight: normal;">(std::vector&lt;double&gt;)</i>

.. |intVecType| raw:: html

   <i style="font-weight: normal;">(std::vector&lt;int&gt;)</i>

The primary interaction is a 2-to-2 collision (projectile + target → ejectile +
residue), where the projectile (target) is the lighter (heavier) initial-state
particle. Four-momenta are evaluated in the laboratory frame (where the target
particle is at rest). The ejectile (residue) is the lighter (heavier)
final-state particle. De-excitation products (γ-rays, neutrons, etc.) are the
particles produced after the primary reaction.

*Projectile*

pdgv |intType|
  PDG code of the projectile.

Ev |doubleType|
  Projectile total energy (MeV).

KEv |doubleType|
  Projectile kinetic energy (MeV).

pxv, pyv, pzv |doubleType|
  Projectile 3-momentum components (MeV).

*Target*

pdgt |intType|
  PDG code of the target.

Mt |doubleType|
  Target mass (MeV).

*Ejectile*

pdgl |intType|
  PDG code of the ejectile.

El |doubleType|
  Ejectile total energy (MeV).

KEl |doubleType|
  Ejectile kinetic energy (MeV).

pxl, pyl, pzl |doubleType|
  Ejectile 3-momentum components (MeV).

*Residue*

pdgr |intType|
  PDG code of the residue.

Er |doubleType|
  Residue total energy (MeV).

KEr |doubleType|
  Residue kinetic energy (MeV).

pxr, pyr, pzr |doubleType|
  Residue 3-momentum components (MeV).

*Residue state immediately following the primary interaction*

Ex |doubleType|
  Excitation energy of the residue immediately after the hard
  interaction and before any de-excitations (MeV).

twoJ |intType|
  Two times the spin J of the residue.

parity |intType|
  Parity of the residue (±1).

*De-excitation products*

np |intType|
  Number of de-excitation products (γ-rays, neutrons, protons, etc.)
  produced during the nuclear de-excitation cascade.

pdgp |intVecType|
  PDG codes of the de-excitation products.

Ep, KEp |doubleVecType|
  Total energies and kinetic energies of the de-excitation products (MeV).

pxp, pyp, pzp |doubleVecType|
  3-momentum components of the de-excitation products (MeV).

tp |doubleVecType|
  Production time for each de-excitation product (seconds). The primary
  interaction occurs at time :math:`t = 0`.

*Cross section and process information*

.. |xsecTreeUnits| raw:: html

   10<sup>&minus;42</sup> cm<sup>2</sup> per atom

xsec |doubleType|
  Flux-averaged total cross section for the run (|xsecTreeUnits|).
  This value is the same for every event in a given run; it is stored
  per-row for convenience.

proc |intType|
  Process ID indicating the kind of primary interaction that occurred in
  the current event. The recognized codes are given in the table below.

  .. list-table::
     :header-rows: 1
     :class: proc-table
     :widths: 8 24 68

     * - ID
       - Name
       - Description
     * - 100
       - vCC-discrete
       - Charged-current neutrino scattering on a nucleus
         that induces a transition to a discrete (bound)
         nuclear energy level
     * - 101
       - vCC-continuum
       - Charged-current neutrino scattering on a nucleus
         that induces a transition to the unbound continuum
         at excitation energies above the particle-emission
         threshold
     * - 110
       - anti-vCC-discrete
       - Same as vCC-discrete (100) but with an incident
         antineutrino
     * - 111
       - anti-vCC-continuum
       - Same as vCC-continuum (101) but with an incident
         antineutrino
     * - 150
       - NC-discrete
       - Neutral-current (anti-)neutrino scattering on a
         nucleus that populates a discrete (bound) energy
         level of the outgoing nucleus
     * - 151
       - NC-continuum
       - Neutral-current (anti-)neutrino scattering on a
         nucleus that populates the unbound continuum
     * - 700
       - v-e
       - Elastic scattering of (anti-)neutrinos on atomic
         electrons
     * - 800
       - standalone-decay
       - Standalone nuclear de-excitation
         (``decay`` command; no primary interaction)

*Event weights*

cv_weight |doubleType|
  The central-value (``CV``) event weight. For normal event generation with
  MARLEY, the central-value weight is always unity.

other_weights |doubleVecType|
  Any additional weights assigned during a ``generate`` or ``reweight`` job.
  These appear in the same order as the strings in the
  ``MARLEY_other_weight_names`` vector mentioned above. The ``other_weights``
  vector may be empty if ``MARLEY_other_weight_names`` is also empty.

Legacy MARLEY v1 format
-----------------------

MARLEY can produce event files in the native text format used in v1.2.1 and
earlier via ``marley convert --output-format legacy``. This is a one-way
conversion provided for backward compatibility; the legacy format is deprecated
and cannot be produced by the ``generate`` command.

HEPEVT
------

The legacy `HEPEVT <https://home.fnal.gov/~mrenna/lutp0613man2/node49.html>`__
output format was officially supported in the MARLEY v1 release series and is
still available as an option for the ``convert`` command. The description
presented here covers only those aspects of the HEPEVT format needed to
interpret the output of MARLEY. Further details are available on pages 327–330
of `this document <https://doi.org/10.5170/CERN-1989-008-V-3>`__.

A HEPEVT-format output file consists of one or more text-based event records.
Each of these records begins with the header

::

  NEVHEP NHEP

where ``NEVHEP`` is the event number (untracked by MARLEY and thus always set
to zero) and ``NHEP`` is the number of particles in the event. The header is
followed by ``NHEP`` lines, each representing a single particle. These have the
format

::

  ISTHEP IDHEP JMOHEP1 JMOHEP2 JDAHEP1 JDAHEP2 PHEP1 PHEP2 PHEP3 PHEP4 PHEP5 VHEP1 VHEP2 VHEP3 VHEP4

where ``ISTHEP`` is an integer code identifying the particle status and
``IDHEP`` is the particle's PDG code. In agreement with the HEPEVT standard,
MARLEY uses status code 1 for the final-state particles and 3 for the
initial-state particles. The ``JMOHEP1``, ``JMOHEP2``, ``JDAHEP1``, and
``JDAHEP2`` entries record the indices (between 1 and ``NHEP``, inclusive) of
particles in the event record that correspond to the first mother, second
mother, first daughter, and last daughter of the current particle,
respectively. These indices are set to zero in cases where they do not apply
(e.g., a particle with no daughters will have ``JDAHEP1`` = ``JDAHEP2`` = 0).
Entries ``PHEP1`` through ``PHEP3`` record the x-, y-, and z-components of the
particle 3-momentum, while ``PHEP4`` gives the total energy and ``PHEP5`` gives
the particle mass (all in GeV). Entries ``VHEP1`` through ``VHEP3`` store the
x, y, and z positions of the particle production vertex (mm), and ``VHEP4``
gives the production time (mm/c).

In addition to the initial- and final-state particles, MARLEY adds a dummy
particle with ``ISTHEP`` = 11 to each HEPEVT event record. All data fields are
zero for this particle except for (1) ``JMOHEP1``, which contains the nuclear
spin multiplied by two, (2) ``JMOHEP2``, which reports the parity of the nucleus
as an integer, (3) ``PHEP4``, which gives the excitation energy of the nucleus
(MeV), and (4) ``PHEP5``, which records the flux-averaged total cross section in
units of MeV:sup:`-2` per atom. The excitation energy, spin, and parity values
in the HEPEVT record refer to the nuclear state that is formed after the primary
interaction but before any de-excitations have occurred.

.. literalinclude:: _static/example.hepevt
   :name: hepevt_format_example

The listing above shows an example MARLEY output file in HEPEVT format.
