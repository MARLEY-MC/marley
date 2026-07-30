=========================
MARLEY command reference
=========================

This page documents each of the subcommands available in the ``marley``
executable. For a quick overview, run ``marley help`` from the command line.
Detailed help for a specific command is available via
``marley help <command>``.

.. contents:: Commands
   :local:
   :depth: 1

--------
generate
--------

Usage::

  marley generate CONFIG_FILE

Generate Monte Carlo neutrino interaction events according to the settings in
``CONFIG_FILE``, a MARLEY job configuration file written in a
JSON-like format.

This is the default command: ``marley CONFIG_FILE`` is equivalent to
``marley generate CONFIG_FILE``.

Examples
~~~~~~~~

Run a simulation using the annotated example configuration file after
sourcing the setup script::

  source setup_marley.sh
  marley generate examples/config/annotated.js

The output file (``events.hepmc3`` by default) will be created in the working
directory.

To override the default output file name, add ``output`` settings to the
configuration file or use ``marley convert`` afterward. See
:doc:`getting_started` for more details on writing configuration files.

-----
print
-----

Usage::

  marley print [FORMAT] INPUT_FILES...

Print one or more MARLEY event files in a human-readable format. This command
reads the native HepMC3 output files and displays their contents on the
terminal. It does not modify the input files.

The optional ``FORMAT`` argument selects the output style:

``pretty`` (default)
  A rich human-readable format that presents event data in an organized table
  with particle labels, energies, and momenta.

``hepmc3``
  Raw HepMC3 text representation streamed directly to stdout. Useful for
  inspecting the underlying HepMC3 event structure.

``legacy``
  A summary format similar to the event display used in MARLEY v1.2.1 and
  earlier. Provided for backward compatibility.

Multiple input files may be specified. They are concatenated in the order
given.

Examples
~~~~~~~~

Print events in the default pretty format::

  marley print events.hepmc3

Show the HepMC3 text representation of two event files::

  marley print --format hepmc3 run1.hepmc3 run2.hepmc3

Use the legacy v1-style summary::

  marley print --format legacy events.hepmc3

--------
reweight
--------

Usage::

  marley reweight CONFIG_FILE INPUT_FILE

Reweight previously generated MARLEY events using weight-calculation settings
from a configuration file. The input event file (in HepMC3 format) is read,
and for each event a weight is computed and written to the output along with
the original event data.

The configuration file should contain a ``weights`` section specifying the
weight calculator(s) to apply, along with any associated parameters. The most
commonly used weight calculator applies systematic variations to tabulated
nuclear matrix elements.

Examples
~~~~~~~~

Reweight an existing event file using a weight configuration::

  marley reweight weight_config.js events.hepmc3

The output file name and format are controlled by the ``output`` settings in
the weight configuration file, following the same conventions as the
``generate`` command.

----
xsec
----

Usage::

  marley xsec -o OUTPUT_FILE CONFIG_FILE

Tabulate the total cross section versus projectile kinetic energy using the
settings in ``CONFIG_FILE``. The results are written to ``OUTPUT_FILE`` as a
two-column text file (energy in MeV, cross section in cm\ :sup:`2` per atom).

The ``-o`` flag is required.

An optional ``xsec`` section may be included in the configuration file
to control the cross-section table parameters. If omitted, defaults are
used for all settings. The supported keys within the ``xsec`` section are:

``KEmin``
  Minimum projectile kinetic energy (MeV). Default: 0.0.

``KEmax``
  Maximum projectile kinetic energy (MeV). Default: 100.0.

``steps``
  Number of equally-spaced energy steps between the minimum and maximum
  kinetic energies. Default: 10000.

``pdg``
  PDG code of the projectile for which to tabulate cross sections.
  Supported values include 12 (ν\ :sub:`e`), 14 (ν\ :sub:`μ`), -12
  (anti-ν\ :sub:`e`), -14 (anti-ν\ :sub:`μ`), etc. Default: 12.

The configuration file must also contain the ``reactions`` and ``source``
keys needed by the Generator (see :doc:`getting_started`). The
``generate`` key is ignored by the ``xsec`` command.

Examples
~~~~~~~~

Generate a cross-section table from a reaction configuration::

  marley xsec -o xsec_table.txt examples/config/minimal.js

Preview the first few lines of the output::

  head -5 xsec_table.txt
  # Projectile kinetic energy (MeV)  Total cross section (cm^2/atom)
  5.00000000000000000e+00            0.00000000000000000e+00
  7.50000000000000000e+00            1.23456789012345679e-40
  1.00000000000000000e+01            9.87654321098765432e-40

-----
decay
-----

Usage::

  marley decay CONFIG_FILE

Simulate stand-alone nuclear de-excitations according to the decay settings
in ``CONFIG_FILE``. This command does not simulate a primary neutrino
interaction; instead, it starts from a user-specified nuclear excited state
and simulates its de-excitation cascade using MARLEY's nuclear de-excitation
models.

This is useful for studying the de-excitation behavior of specific nuclear
levels in isolation, or for generating de-excitation-only event samples.

The configuration file uses the same JSON-like format as the ``generate``
command, but should specify the initial nuclear state (excitation energy,
spin, parity) rather than a reaction input file.

Examples
~~~~~~~~

Run a stand-alone decay simulation::

  marley decay my_decay_config.js

---------
summarize
---------

Usage::

  marley summarize -o OUTPUT_FILE INPUT_FILE...

Convert MARLEY event files into a ROOT TTree summary file suitable for
analysis with ROOT C++ macros or Python (PyROOT). This produces a
flat ntuple where each event is a single row in the tree, with branches for
the projectile, target, ejectile, residue, and de-excitation products.

The ``-o OUTPUT_FILE`` flag is required. Multiple input files may be
specified; they will be concatenated in the order given.

This command requires a ROOT-enabled build of MARLEY. See
:doc:`getting_started` for build instructions.

The structure of the output TTree is described in the :doc:`interpret_output`
documentation.

Examples
~~~~~~~~

Summarize a generated event file into a flat ROOT ntuple::

  marley summarize -o summary.root events.hepmc3

Combine two event files into a single summary tree::

  marley summarize -o combined.root run1.hepmc3 run2.hepmc3

-------
convert
-------

Usage::

  marley convert [--output-format FORMAT] -o OUTPUT_FILE INPUT_FILES...

Convert MARLEY event files between supported output formats. This command
reads one or more input files in the native HepMC3 format and writes them in
the requested output format.

Options:

``-o OUTPUT_FILE``
  Required: path to the output file.

``--output-format FORMAT``
  Output format selector. Supported values:

  ``ascii`` (default)
    HepMC3 standard ASCII format (same as the default ``generate`` output).

  ``root``
    HepMC3-based ROOT format using ``GenEventData`` structures.
    Requires a ROOT-enabled build of MARLEY.

  ``legacy``
    MARLEY v1.2.1 native text format (one-way conversion).

  ``hepevt``
    Legacy HEPEVT text format (one-way conversion).

``-f, --force``
  Overwrite the output file without prompting.

The ``legacy`` and ``hepevt`` formats are one-way conversions from the
current HepMC3-based format to the deprecated output formats used in MARLEY
v1.2.1 and earlier.

Examples
~~~~~~~~

Convert an event file to the legacy v1 text format::

  marley convert --output-format legacy -o events_legacy.txt events.hepmc3

Convert to HEPEVT format::

  marley convert --output-format hepevt -o events.hepevt events.hepmc3

Convert multiple files to ROOT format (requires ROOT build)::

  marley convert --output-format root -o combined.root run1.hepmc3 run2.hepmc3

Overwrite an existing output file without a confirmation prompt::

  marley convert --output-format ascii -o output.hepmc3 input.hepmc3 -f

----
help
----

Usage::

  marley help [COMMAND]

Display the top-level help message (listing all commands) or detailed help
for a specific command. Equivalent to ``marley --help``
or ``marley <command> --help``.

Examples
~~~~~~~~

Show the top-level help::

  marley help

Show detailed help for the ``generate`` command::

  marley help generate

-------
version
-------

Usage::

  marley version

Print version information, including the MARLEY release version, the git
revision hash, and the build date. Equivalent to ``marley --version``.

Example::

  $ marley version
  MARLEY v2.0.0 (git revision abc1234, built Jul 19 2026)
