=======================
Interpreting the output
=======================

This page provides a guide to the contents of the output files produced by the
``marley`` executable. Following a brief description of the *PDG codes* used to
identify particle types in MARLEY, documentation for each available output
format is given below.

PDG codes
^^^^^^^^^

The `Particle Data Group <https://pdg.lbl.gov/index.html>`__ (PDG) has defined a
standard numbering scheme for representing particle species in Monte Carlo
event generators. Each kind of particle is assigned a unique positive integer
as an identifier. The corresponding antiparticle is assigned a negative integer
with the same absolute value. A full description of the numbering scheme is
available `here <https://pdg.lbl.gov/current/mc-particle-id>`__. 

Like nearly all modern particle physics generators, MARLEY adopts the integer
*PDG codes* for particle identification and uses them both internally and in
output files. For convenience, a table of the PDG codes most relevant for
MARLEY is given below.

.. raw:: html
   
     <head>
     <style>
       table.mytable {
         border-collapse: separate;
         border-spacing: 0 5px;
       }
       table.mytable th, td {
         text-align: center;
       }
     </style>
     </head>

   <table align="center" style="width:30%" class="mytable">

     <tr>
       <th>PDG code</th>
       <th>Particle</th>
     </tr>

     <tr>
       <td>11</td>
       <td>e<sup>&minus;</sup></td>
     </tr>

     <tr>
       <td>12</td>
       <td>&nu;<sub>e</sub></td>
     </tr>

     <tr>
       <td>13</td>
       <td>&mu;<sup>&minus;</sup></td>
     </tr>

     <tr>
       <td>14</td>
       <td>&nu;<sub>&mu;</sub></td>
     </tr>

     <tr>
       <td>15</td>
       <td>&tau;<sup>&minus;</sup></td>
     </tr>

     <tr>
       <td>16</td>
       <td>&nu;<sub>&tau;</sub></td>
     </tr>

     <tr>
       <td>22</td>
       <td>&gamma;</td>
     </tr>

     <tr>
       <td>2112</td>
       <td>n</td>
     </tr>

     <tr>
       <td>2212</td>
       <td>p</td>
     </tr>

     <tr>
       <td>1000010020</td>
       <td>d</td>
     </tr>

     <tr>
       <td>1000010030</td>
       <td>t</td>
     </tr>

     <tr>
       <td>1000020030</td>
       <td>h</td>
     </tr>


     <tr>
       <td>1000020040</td>
       <td>&alpha;</td>
     </tr>

  </table>

In general, a PDG code of the form 100ZZZAAA0 represents a nuclide with proton
number Z and mass number A. For example, :superscript:`40`\Ar is represented by
the PDG code 1000180400.

Output file formats
^^^^^^^^^^^^^^^^^^^

The ``marley generate`` and ``marley decay`` commands produce event files in
the HepMC3 standard `ASCII <https://hepmc.org>`__ format by default. The
``marley convert`` command can translate between this format and several
additional formats for backward compatibility. Descriptions of each available
format are given below.

HepMC3 (default)
----------------

The native output format for MARLEY v2.0.0 and later is the standard ASCII
representation defined by the `HepMC3 <https://hepmc.org>`__ event record
library (``WriterAscii`` class). This format is produced by default by the
``marley generate`` and ``marley decay`` commands, and by ``marley convert``
when no ``--output-format`` argument is supplied.

Each HepMC3 output file contains:

* A run-info block with metadata attributes following the
  `NuHepMC <https://github.com/NuHepMC/standard>`__ conventions,
  including generator name, version, and citation identifiers.
* One or more ``GenEvent`` records, each representing a single neutrino
  interaction (or nuclear de-excitation). Events are structured as a
  graph of vertices and particles, with each vertex recording the
  production and end points of particles.

Within each event, the primary interaction is modeled as a 2-to-2
scattering::

  projectile + target → ejectile + residue

where the projectile is the lighter of the two initial-state particles,
the target is the heavier (and is at rest in the lab frame), the ejectile
is the lighter final-state particle, and the residue is the heavier.
De-excitation products (γ-rays, neutrons, etc.) are stored as daughters
of the final-state vertices.

The full specification of the HepMC3 format is maintained by the HepMC3
project (DOI `10.1016/j.cpc.2020.107310
<https://doi.org/10.1016/j.cpc.2020.107310>`__). Users who wish to
process MARLEY output programmatically are encouraged to use the HepMC3
library's API rather than parsing the text representation directly.

Legacy MARLEY v1 format
-----------------------

MARLEY can produce event files in the native text format used in
v1.2.1 and earlier via ``marley convert --output-format legacy``.
This is a one-way conversion provided for backward compatibility;
the legacy format is no longer the default and will not be produced
by ``marley generate``. Users who need to read legacy-format files
should use ``marley convert`` or ``marley print`` to view them.


HEPEVT
------

The legacy `HEPEVT <https://home.fnal.gov/~mrenna/lutp0613man2/node49.html>`__
format is designed for interfacing event generators with each other and with
other software. Compatibility with this event format is maintained in many
modern high energy physics software libraries. The description presented here
covers only those aspects of the HEPEVT format needed to interpret the output
of MARLEY. Further details are available on pages 327–330 of `this document
<https://doi.org/10.5170/CERN-1989-008-V-3>`__.

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

Because MARLEY currently treats all nuclear de-excitations as instantaneous and
does not perform any particle tracking, ``VHEP1`` through ``VHEP4`` are always
identically zero in HEPEVT output files. Intermediate de-excitation steps are
also not currently stored in the event record, so ``JMOHEP1``, ``JMOHEP2``,
``JDAHEP1``, and ``JDAHEP2`` are also identically zero in most cases.

In addition to the initial- and final-state particles, MARLEY adds a dummy
particle with ``ISTHEP`` = 11 to each HEPEVT event record. All data fields are
zero for this particle except for (1) ``JMOHEP1``, which contains the nuclear
spin multiplied by two, (2) ``JMOHEP2``, which reports the parity of the
nucleus as an integer, (3) ``PHEP4``, which gives the excitation energy of the
nucleus (MeV), and (4) ``PHEP5``, which records the flux-averaged total cross
section in units of |InverseMeVSquared| per atom. As is the case for the
HepMC3 format, the excitation energy, spin, and parity values refer to the
nuclear state that is formed after the primary scattering reaction but before
any de-excitations have occurred.

.. literalinclude:: _static/example.hepevt
   :name: hepevt_format_example

The listing above shows an example MARLEY output file in HEPEVT format. The
same two events from the HEPEVT-format example file
are used for easy comparison of the formats.

ROOT
----

If MARLEY has been built with ROOT support (see the :doc:`getting_started`
page for build instructions), events can be stored in ROOT's compressed
binary format using the ``root`` output format of ``marley generate`` or
``marley convert --output-format root``.

MARLEY v2.0.0 uses the HepMC3 ``WriterRoot`` class to produce ROOT output.
The resulting files contain a ``TTree`` named ``hepmc3_event_tree`` with a
single branch of ``GenEventData`` objects. This format does not require
MARLEY-specific class dictionaries; the standard HepMC3 ROOT dictionary is
used instead.

The ``mroot`` helper script from MARLEY v1.2.1 and earlier is no longer
provided. To analyze ROOT output files, load the HepMC3 ROOT dictionary:

* From a ROOT C++ macro or the ROOT prompt::
  
    R__LOAD_LIBRARY(libHepMC3rootIO)
  
* From PyROOT (requires ``pip install hepmc3`` or equivalent)::
  
    import ROOT
    ROOT.gSystem.Load("libHepMC3rootIO")

Example analysis macros are available in the ``examples/macros/`` folder
of the MARLEY source distribution.

Conventions used in the event objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When analyzing ``marley::Event`` objects (or their alternative "flat"
representation discussed `below <#flat-root-files>`__), it is helpful to be
aware of MARLEY's nomenclature for the particles involved in each event's 2 → 2
primary interaction
|genericReaction|

Here the particles 𝑎, 𝑏, 𝑐, and 𝑑 are labeled as, respectively, the
*projectile*, *target*, *ejectile*, and *residue*. Where a mass difference
exists, MARLEY chooses the projectile (ejectile) to be the lighter of the two
particles in the initial (final) state. Otherwise, the choice is arbitrary. All
four-vector components stored in a MARLEY event record are expressed in the
laboratory frame, i.e., the rest frame of the target. If simulation of nuclear
de-excitations is enabled (as it is by default) and the residue is a nucleus,
its 4-momentum and net charge are stored in the event record after it has
reached the ground state.

Metadata
~~~~~~~~

Four pieces of metadata are saved in a ROOT-format output file alongside the
events themselves:

MARLEY_config
  A JSON-format string which stores the contents (except for comments) of the
  job configuration file used to generate the events

MARLEY_state
  A string giving the serialized internal state (obtained using the stream
  insertion operator ``<<``) of the `std::mt19937_64
  <https://en.cppreference.com/w/cpp/numeric/random/mt19937_64>`__ object used to
  generate random numbers when creating the events.

MARLEY_seed
  A string representation of the integer random number seed used to initialize
  the simulation

MARLEY_flux_avg_xsec 
  A `TParameter\<double\>
  <https://root.cern/root/html528/TParameter_double_.html>`__ which stores the
  flux-averaged total cross section (|InverseMeVSquared| per atom)

"Flat" ROOT files
~~~~~~~~~~~~~~~~~

An alternative "flat" form of the ROOT output format is also available which
may be analyzed without the need for the MARLEY class dictionaries. An output
file containing MARLEY events in any of the four standard formats may be
converted into a "flat" ROOT file using the ``marley summarize`` command.
After sourcing the `setup_marley.sh
<getting_started.html#setting-up-the-runtime-environment>`__ script, one may
convert the MARLEY output file ``OLD_EVENTS_FILE`` into a new "flat" ROOT file,
``new_flat_file.root``, via the command

::

  marley summarize -o new_flat_file.root OLD_EVENTS_FILE

.. |doubleType| raw:: html

   <i style="font-weight: normal;">(double)</i>

.. |intType| raw:: html

   <i style="font-weight: normal;">(int)</i>

.. |doubleArrayType| raw:: html

   <i style="font-weight: normal;">(double[np])</i>

.. |intArrayType| raw:: html

   <i style="font-weight: normal;">(int[np])</i>

.. |sumTree| raw:: html

   <u>M</u>ARLEY <u>s</u>ummary <u>t</u>ree


The new file will contain
a single ROOT ``TTree`` called ``mst``
(for |sumTree|) with the following branches:

.. Projectile

pdgv |intType|
  Projectile PDG code

Ev |doubleType|
  Projectile total energy (MeV)

KEv |doubleType|
  Projectile kinetic energy (MeV)

pxv |doubleType|
  Projectile 3-momentum x-component (MeV)

pyv |doubleType|
  Projectile 3-momentum y-component (MeV)

pzv |doubleType|
  Projectile 3-momentum z-component (MeV)

.. Target

pdgt |intType|
  Target PDG code

Mt |doubleType|
  Target mass (MeV)

.. Ejectile

pdgl |intType|
  Ejectile PDG code

El |doubleType|
  Ejectile total energy (MeV)

KEl |doubleType|
  Ejectile kinetic energy (MeV)

pxl |doubleType|
  Ejectile 3-momentum x-component (MeV)

pyl |doubleType|
  Ejectile 3-momentum y-component (MeV)

pzl |doubleType|
  Ejectile 3-momentum z-component (MeV)

.. Residue

pdgr |intType|
  Residue PDG code

Er |doubleType|
  Residue total energy (MeV)

KEr |doubleType|
  Residue kinetic energy (MeV)

pxr |doubleType|
  Residue 3-momentum x-component (MeV)

pyr |doubleType|
  Residue 3-momentum y-component (MeV)

pzr |doubleType|
  Residue 3-momentum z-component (MeV)

.. Residue state before de-excitations

Ex |doubleType|
  Initial residue excitation energy (MeV)

twoJ |intType|
  Two times the initial residue spin

parity |intType|
  Initial residue parity

.. De-excitation products

np |intType|
  Number of de-excitation products

pdgp |intArrayType|
  De-excitation product PDG codes

Ep |doubleArrayType|
  De-excitation product total energies (MeV)

KEp |doubleArrayType|
  De-excitation product kinetic energies (MeV)

pxp |doubleArrayType|
  De-excitation product 3-momentum x-components (MeV)

pyp |doubleArrayType|
  De-excitation product 3-momentum y-components (MeV)

pzp |doubleArrayType|
  De-excitation product 3-momentum z-components (MeV)

.. Flux-averaged total cross section

.. |xsecTreeUnits| raw:: html

   10<sup>&minus;42</sup> cm<sup>2</sup> per atom

xsec |doubleType|
  Flux-averaged total cross section (|xsecTreeUnits|)

.. |InverseMeVSquared| raw:: html

   MeV<sup> &minus;2</sup>

.. |genericReaction| raw:: html

   <p style="text-align: center;"> 𝑎 + 𝑏 → 𝑐 + 𝑑 .</p>
