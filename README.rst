MARLEY (Model of Argon Reaction Low Energy Yields)
==================================================

|platform| |License: GPL v3| |DOI|

|Build Status| |rel| |commits since|

Introduction
------------

.. overview-start

**MARLEY** (Model of Argon Reaction Low Energy Yields) is a Monte Carlo event
generator for neutrino interactions at energies of tens-of-MeV and below. The
current version computes inclusive neutrino-nucleus cross sections using a
hybrid model; a Hartree-Fock Continuum Random Phase Approximation approach is
used for induced transitions to unbound nuclear energy levels, while a
data-driven recipe is applied for transitions to bound states. De-excitations of
the final-state nucleus emerging from the primary interaction are simulated
using a combination of tabulated γ-ray decay schemes and an original
implementation of the Hauser-Feshbach statistical model.

Input files are provided with the code that are suitable for simulating the
charged-current process

.. overview-math-start

|ve40ArCC|

.. overview-math-end

coherent elastic neutrino-nucleus scattering (CEvNS) on spin-zero target
nuclei, and neutrino-electron elastic scattering on any atomic target.
Inclusion of additional reactions and targets is planned for the future.

.. overview-end

MARLEY follows an open-source development model and welcomes contributions of
new input files and code improvements from the community. A partial list of
potential projects for future MARLEY development is available on the developer
documentation `webpage
<https://www.marleygen.org/dev_docs.html#development-wish-list>`__.

Copyright and License
---------------------

.. copyright-start-1

Copyright © 2016-2026 Steven Gardiner gardiner@fnal.gov

MARLEY is distributed under the terms of version 3 of the `GNU General Public
License <https://www.gnu.org/licenses/gpl-3.0-standalone.html>`__ ("GPLv3") as
published by the Free Software Foundation. For the full text of that license,
please see the `COPYING <COPYING>`__ file.

.. copyright-start-2
As a matter of professional courtesy, MARLEY users are also requested to follow
the `MCnet Guidelines <https://www.montecarlonet.org/GUIDELINES>`__ for Event
Generator Authors and Users. Nevertheless, these guidelines are not legally
binding and do not limit your rights guaranteed under the GPLv3.
See the `GUIDELINES <GUIDELINES>`__ file for more details.

Neither the United States nor the United States Department of Energy, nor any of
their employees, makes any warranty, express or implied, or assumes any legal
liability or responsibility for the accuracy, completeness, or usefulness of any
data, apparatus, product, or process disclosed, or represents that its use would
not infringe privately owned rights.

Citing MARLEY
-------------

.. citing-start

If you refer to MARLEY in academic work, please **always cite** the following
reference:

S. Gardiner, Simulating low-energy neutrino interactions with MARLEY,
`Comput. Phys. Commun. 269, 108123
<https://doi.org/10.1016/j.cpc.2021.108123>`__,
`arXiv:2101.11867 [nucl-th] <https://arxiv.org/abs/2101.11867>`__ (2021).

.. citing-math-start

In publications which use the recommended physics configuration for
charged-current νₑ-⁴⁰Ar scattering (see the REACTION INPUT FILES section of
examples/config/annotated.js), please also cite the paper describing the MARLEY
v2 model of this reaction:

S. Gardiner *et al.*, Continuum contribution to charged-current absorption of
low-energy νₑ on ⁴⁰Ar, `arXiv:2604.26801 [hep-ph]
<https://arxiv.org/abs/2604.26801>`__ (2026).

.. citing-math-end

Providing a citation for the MARLEY code itself is also encouraged and
appreciated. To maximize reproducibility of published calculations, such
citations should include the digital object identifier (DOI) associated with the
code release that was used. The DOIs for recent versions of MARLEY are listed on
the GitHub `releases webpage <https://github.com/MARLEY-MC/marley/releases>`__
and in the right-hand column of the Zenodo "concept DOI" `webpage
<https://doi.org/10.5281/zenodo.3901933>`__.

For convenience, recommended `BibTeX <https://www.bibtex.org/>`__ citations to
use for the latest MARLEY release are given in the `CITATION.bib
<CITATION.bib>`__ file.

.. citing-end

If you use the default nuclear structure data files (strongly recommended) for
published calculations, please also give proper attribution to the developers
of the TALYS nuclear code (see `data/structure/README.md
<data/structure/README.md>`__ for more information).

Getting Started
---------------

.. getting-started-start1

MARLEY is regularly tested on both Linux and macOS platforms and is expected to
work in any Unix-like environment in which the prerequisites are installed.
Building and running MARLEY on Windows is not currently supported.

Prerequisites
~~~~~~~~~~~~~

There are two prerequisites needed to build MARLEY:

.. getting-started-end1

.. class:: open

.. getting-started-start2

*  A C++17-compliant compiler. The following compilers are officially
   supported:

   -  `GNU Compiler Collection <https://gcc.gnu.org>`__ (GCC) ≥ 9.1.0

   -  `Clang <https://clang.llvm.org>`__ ≥ 9.0.0

*  `GNU Make <https://www.gnu.org/software/make/>`__

On both Linux and macOS, these prerequisites will likely be available through
the standard package manager. Note that when `CMake <https://cmake.org/>`__ is
available on the host system, it is used by default. However, an equivalent
build configuration using only GNU Make is also provided for user convenience.

Although it is not required in order to build or use MARLEY, the popular `ROOT
<https://root.cern.ch>`__ data analysis framework provides convenient tools for
plotting and analyzing simulation results. Users who wish to use the optional
interface between the two codes should ensure that ROOT is installed before
building MARLEY. At build time, the optional MARLEY interface to ROOT is
enabled automatically if the ``root-config`` script is present on the system
``PATH``.

MARLEY has two additional dependencies that are both optional: the `GNU
Scientific Library <https://www.gnu.org/software/gsl/>`__ (GSL) and the `HepMC3
<https://gitlab.cern.ch/hepmc/HepMC3>`__ event record library. If these
dependencies are not detected on the host system at build time, then built-in
versions will be used instead.

.. getting-started-end2

.. getting-started-start3

Building MARLEY
~~~~~~~~~~~~~~~

To build the code, run ``make`` from the top-level MARLEY directory:

::

    make

The top-level ``Makefile`` will auto-detect `CMake <https://cmake.org/>`__ and
use it as the build backend if available; otherwise it falls back to the
included `GNU Make <https://www.gnu.org/software/make/>`__ recipe
(``make/build.mk``). The ``build/`` directory is created automatically by
either backend and is removed by running ``make clean``.

If the build is successful, then the ``marley`` executable will be located at
``build/bin/marley``. Running it without arguments

::

    build/bin/marley

should produce the following output:

::

    Usage: marley <command> [options]

    Commands:
      convert     Convert event files between supported formats
      decay       Simulate nuclear de-excitations
      generate    Generate Monte Carlo events
      help        Show this help message or help for a specific command
      print       Print existing events in a human-readable format
      reweight    Reweight previously generated events
      summarize   Create a ROOT TTree summary of event files [requires ROOT]
      version     Print version information
      xsec        Tabulate total cross section vs. projectile kinetic energy

    Options:
      -h, --help     Show top-level help
      -v, --version  Print version information

    Run 'marley help <command>' or 'marley <command> --help' for details.
    MARLEY home page: <https://www.marleygen.org>

From the top-level Makefile, the user can optionally direct the build system
to ignore CMake (thus falling back to a pure GNU Make recipe) as well as
any of the optional dependencies. Invoking ``make`` with the settings

::

    make IGNORE_CMAKE=1 IGNORE_ROOT=1 IGNORE_GSL=1 IGNORE_HEPMC3=1

will bypass CMake when building the code, disable the ROOT interface even if a
ROOT installation is successfully detected, and force the use of built-in
versions of GSL and HepMC3 even if system installations are available for both
of these libraries. Any combination of these ``make`` options may be used
in any order according to the user's preferences.

Setting up the runtime environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``marley`` executable relies on the system environment variable ``MARLEY``
to store the full path to the root folder of the source code. This variable may
be set automatically by sourcing the ``setup_marley.sh`` Bash script:

::

    source setup_marley.sh

For user convenience, this script also adds ``build/bin`` to the system
``PATH`` and adds ``build/lib`` to ``LD_LIBRARY_PATH`` and, on macOS,
``DYLD_LIBRARY_PATH``. After sourcing the setup script, the ``marley``
command may be run from any directory.

If generation of events is attempted without setting the ``MARLEY`` environment
variable first, then MARLEY will halt after printing the error message

::

    [ERROR]: The MARLEY environment variable is not set. Please set it (e.g.,
    by sourcing the setup_marley.sh script) and try again.

Generating events
~~~~~~~~~~~~~~~~~

The ``marley`` executable allows the user to adjust simulation parameters
via job configuration files written in a `JSON
<https://www.json.org/json-en.html>`__-like format. The name of the
configuration file to use appears as the first argument after the
``generate`` command:

::

  marley generate CONFIG_FILE

To generate events using an example configuration file, execute the following
command after sourcing the ``setup_marley.sh`` script:

::

    marley generate examples/config/annotated.js

The program will display the MARLEY logo and diagnostic messages as it runs the
simulation. When the program terminates, a new file named ``events.hepmc3`` will
be present in the working directory. This file contains the generated events in
the standard ASCII representation of the `HepMC3
<https://doi.org/10.1016/j.cpc.2020.107310>`__ data format.

The ``annotated.js`` configuration file mentioned above is heavily commented
with explanations of the most commonly-used input parameters. Reading it serves
as a good next step for new users. When you are ready to start writing your own
configuration files, editing a copy of ``examples/config/COPY_ME.js`` is
recommended.

.. getting-started-end3

Core Developers
---------------

.. class:: open

- **Steven Gardiner** - `sjgardiner <https://github.com/sjgardiner>`__

See also the list of `contributors
<https://github.com/MARLEY-MC/marley/contributors>`__ who participated in this
project.

Website
-------

Further documentation for the latest version of MARLEY may be found on the
official website at `https://www.marleygen.org/ <https://www.marleygen.org>`__.

Acknowledgements
----------------

Special thanks go to

.. class:: open

- The `TALYS <https://talys.eu>`__ developers (Arjan Koning, Stéphane
  Hilaire, and Marieke Duijvestijn) for sharing their nuclear structure data

- Zero Anixter for providing an illustration of Bob Marley to be used
  in the official MARLEY logo

.. |ve40ArCC| raw:: html

   <p align="center">&nu;<sub>e</sub>&nbsp;+&nbsp;<sup>40</sup>Ar&nbsp;&rarr;
   &nbsp;e<sup>&minus;</sup>&nbsp;+&nbsp;<sup>40</sup>K<sup>&ast;</sup>,</p>

.. |platform| image:: https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey

.. |License: GPL v3| image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

.. |DOI| image:: https://img.shields.io/badge/DOI-10.5281%2Fzenodo.3901933-blue
   :target: https://doi.org/10.5281/zenodo.3901933

.. |Build Status| image:: https://github.com/MARLEY-MC/marley/actions/workflows/ci.yml/badge.svg?branch=main
   :target: https://github.com/MARLEY-MC/marley/actions/workflows/ci.yml

.. |rel| image:: https://img.shields.io/github/v/release/MARLEY-MC/marley?include_prereleases
   :target: https://github.com/MARLEY-MC/marley/releases

.. |commits since| image:: https://img.shields.io/github/commits-since/MARLEY-MC/marley/latest/main
   :target: https://github.com/MARLEY-MC/marley/commits/main
