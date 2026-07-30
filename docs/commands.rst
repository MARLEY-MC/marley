=================
Command reference
=================

This page documents each of the commands recognized by the ``marley``
executable. For a quick overview, run ``marley help`` from the command line.
Detailed help for a specific command is available via ``marley help <command>``
or ``marley <command> --help``.

.. raw:: html

   <div style="margin-bottom: 2rem;"></div>

.. contents:: Commands
   :local:
   :depth: 1

-------
convert
-------

Usage::

  marley convert [--output-format FORMAT] -o OUTPUT_FILE INPUT_FILES...

This command reads one or more input files containing HepMC3-format MARLEY
events and writes them in the requested output format.

Options:

``-o OUTPUT_FILE``
  Required: path to the output file.

``--output-format FORMAT``
  Output format selector. Supported values:

  ``ascii``
    HepMC3 standard ASCII format (same as the default ``generate`` output).

  ``root``
    HepMC3 standard `ROOT <https://root.cern>`__ format. Requires a
    ROOT-enabled build of MARLEY.

  ``legacy``
    MARLEY v1 native text format (one-way conversion).

  ``hepevt``
    Legacy HEPEVT text format (one-way conversion).

If ``--output-format`` is omitted, the format is inferred from the
output filename: ``.root`` → ``root``, otherwise ``ascii``.

``-f, --force``
  Overwrite the output file without prompting to confirm.

The ``legacy`` and ``hepevt`` formats are one-way conversions from the current
HepMC3-based format to the deprecated output formats used in MARLEY v1.2.1 and
earlier.

Examples
~~~~~~~~

Convert an event file to the legacy MARLEY v1 text format::

  marley convert --output-format legacy -o events_legacy.txt events.hepmc3

Convert to HEPEVT format::

  marley convert --output-format hepevt -o events.hepevt events.hepmc3

Convert multiple files to ROOT format (requires building MARLEY against ROOT)::

  marley convert --output-format root -o combined.root run1.hepmc3 run2.hepmc3

Overwrite an existing output file without a confirmation prompt::

  marley convert -f --output-format ascii -o output.hepmc3 input.hepmc3

-----
decay
-----

Usage::

  marley decay CONFIG_FILE

Simulate stand-alone nuclear de-excitations according to the settings in
``CONFIG_FILE``. This command does not simulate a primary neutrino interaction;
instead, it starts from a user-specified nuclear excited state and simulates its
de-excitation cascade. The treatment of de-excitation physics is the same as for
neutrino scattering events.

The example ``marley decay`` configuration file
(``examples/config/decay_example.js``) provides documentation for all relevant
settings, some of which are shared with the ``generate`` command.

Example
~~~~~~~

Run a simulation of nuclear de-excitations using the example configuration::

  marley decay examples/config/decay_example.js

--------
generate
--------

Usage::

  marley generate CONFIG_FILE

Generate neutrino interaction events according to the settings in
``CONFIG_FILE``, a MARLEY job configuration file written in a JSON-like format.
The annotated example job configuration file (``examples/config/annotated.js``)
provides detailed documentation for settings that may be used with the
``generate`` command.

This is the default command: ``marley CONFIG_FILE`` is equivalent to
``marley generate CONFIG_FILE``.

Example
~~~~~~~

Run a simulation using the annotated example configuration file after
sourcing the setup script::

  source setup_marley.sh
  marley generate examples/config/annotated.js

The output file (``events.hepmc3`` by default) will be created in the working
directory. To override the default output file name, edit the ``output``
settings in the job configuration file.

----
help
----

Usage::

  marley help [COMMAND]

Display the top-level help message (listing all commands) or detailed help
for a specific command. Equivalent to ``marley --help``
or ``marley <command> --help``.

Example
~~~~~~~

Show detailed help for the ``generate`` command::

  marley help generate

-----
print
-----

Usage::

  marley print [FORMAT] INPUT_FILES...

Print events from one or more MARLEY output files to the terminal in a
human-readable format. The contents of the files are not altered in any way.

The optional ``FORMAT`` argument selects the printing style:

``pretty`` (default)
  The recommended format for visual inspection of MARLEY events. A summary
  of the event is shown first, followed by detailed information
  about each interaction vertex. The printout concludes with a table of the
  final-state particles, their 4-momenta, and cross-section values for the
  sampled energy of the projectile.

``hepmc3``
  Raw text representation of the `HepMC3
  <https://doi.org/10.1016/j.cpc.2020.107310>`__ event format. This is obtained
  using a C++ stream operator on the ``HepMC3::GenEvent`` object defined in the
  official `HepMC3 library <https://gitlab.cern.ch/hepmc/HepMC3>`__.

``legacy``
  The event summary format produced by the ``marprint`` executable distributed
  with MARLEY v1.2.1 and earlier. This format is kept for backward
  compatibility.

Multiple input files may be specified. Events are printed from each file in the
order that they appear.

Examples
~~~~~~~~

Print events in the default "pretty" format::

  marley print events.hepmc3

Show the HepMC3 text representation of two event files::

  marley print hepmc3 run1.hepmc3 run2.hepmc3

Use the legacy summary format from MARLEY v1::

  marley print legacy events.hepmc3

--------
reweight
--------

Usage::

  marley reweight CONFIG_FILE INPUT_FILE...

Compute new weights for an existing sample of MARLEY events. One or more
HepMC3-format files containing the original events are provided as input. For
each event, the new weight(s) are computed and appended to the existing list.
The full event records including the new weights are written to the output.

Event reweighting is a common technique for evaluation of interaction model
uncertainties in neutrino experiments. A description of early work on this
subject for the `GENIE <https://genie-mc.org>`__ neutrino event generator
is available `here <https://www.actaphys.uj.edu.pl/R/40/9/2613/pdf>`__.
The MARLEY implementation of event reweighting follows similar principles.

Documentation for the settings used by ``marley reweight`` is available
in an example configuration file (``examples/config/reweight_example.js``).

Example
~~~~~~~

Reweight an existing sample of events using the example configuration file::

  marley reweight examples/config/reweight_example.js events.hepmc3

The output file name(s) and format(s) are controlled by the ``output`` settings
within the ``reweight`` section of the configuration file. The syntax for
configuring the output is the same as for the ``generate`` command.

Notes
~~~~~

* When multiple input files are specified, MARLEY checks the run information
  in each file to verify compatibility. If a discrepancy in the physics
  configuration is found between the files, then the reweighting job will
  terminate with an error message.

* The ``marley reweight`` command can be run multiple times to append new
  weights as many times as desired. However, weight names must be unique across
  all such reweighting jobs. The executable will halt with an error if a new
  weight is requested whose name matches the name of an existing weight.

---------
summarize
---------

Usage::

  marley summarize -o OUTPUT_FILE INPUT_FILE...

Convert one or more MARLEY event files into a `ROOT <https:://root.cern>`__
summary ntuple file. This summary file provides a convenient representation of
the events for subsequent analysis with ROOT using either C++ or `Python
<https://root.cern/manual/python/>`__. The summary ntuple file format is
described :ref:`here <summary-ntuple-format>`.

The ``-o OUTPUT_FILE`` flag is required. Multiple input files may be specified;
they will be concatenated in the order given. If the physics configuration is
found to be inconsistent between the input files, then the executable will halt
with an error message.

This command requires a ROOT-enabled build of MARLEY. See
:doc:`getting_started` for build instructions.

Example
~~~~~~~

Create a summary ROOT ntuple from a HepMC3-format file containing MARLEY events::

  marley summarize -o summary.root events.hepmc3

-------
version
-------

Usage::

  marley version

Print version information for the current build of MARLEY. Equivalent to
``marley --version``.

Example::

  $ marley version
  MARLEY (Model of Argon Reaction Low Energy Yields) 2.0.0
  Copyright (C) 2016-2026 Steven Gardiner
  License: GNU GPL version 3 <http://opensource.org/licenses/GPL-3.0>
  This is free software: you are free to change and redistribute it.

----
xsec
----

Usage::

  marley xsec -o OUTPUT_FILE CONFIG_FILE

Tabulate the total cross section versus projectile kinetic energy using the
settings in ``CONFIG_FILE``. The results are written to ``OUTPUT_FILE`` as a
two-column text file (energy in MeV, cross section in
:math:`10^{-42}\,\mathrm{cm}^{2} / \mathrm{atom}`).

The ``-o`` flag is required.

An optional ``xsec`` section may be included in the configuration file to
control production of the cross-section table. If omitted, defaults are used
for all settings. There are four supported keys within the ``xsec`` section:

``KEmin``
  Minimum projectile kinetic energy (MeV). Default: 0.0.

``KEmax``
  Maximum projectile kinetic energy (MeV). Default: 100.0.

``steps``
  Number of rows in the tabulated grid, spanning :math:`[KE_{min}, KE_{max}]`
  inclusive. Default: 10000.

``pdg``
  :ref:`PDG code <pdg-codes>` of the projectile for which to tabulate cross
  sections. Default: 12 (:math:`\nu_e`).

Example
~~~~~~~~

Generate a cross-section table from the example ``marley generate``
configuration::

  marley xsec -o xsec_table.txt examples/config/annotated.js

The output is a plain two-column text file. Each row gives one kinetic-energy
point (MeV) and the corresponding inclusive total cross section
(:math:`10^{-42}\,\mathrm{cm}^{2} / \mathrm{atom}`). For example, the first few
rows produced using the example command given above are
::

  0 0
  0.010001 4.95196e-05
  0.020002 0.000192741
  0.030003 0.000422255
  0.040004 0.000731378

The tabulated cross section is the inclusive sum over all active reactions in
the configuration. Since the charged-current reaction on
:math:`^{40}\mathrm{Ar}` has a few-MeV threshold, the nonzero cross section
values seen above come from neutrino-electron elastic scattering
(``data/react/ES.react``).
