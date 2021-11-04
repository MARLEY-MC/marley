MARLEY (Model of Argon Reaction Low Energy Yields)
==================================================

|platform| |License: GPL v3| |DOI|

|Build Status| |rel| |commits since|

Introduction
------------

.. overview-start

.. |gamma| unicode:: 0x3B3 .. lowercase gamma

**MARLEY** (Model of Argon Reaction Low Energy Yields) is a Monte Carlo event
generator for neutrino-nucleus interactions at energies of tens-of-MeV and
below. The current version computes inclusive neutrino-nucleus cross sections
employing the *allowed approximation*: the nuclear matrix elements are
evaluated while neglecting Fermi motion and applying the long-wavelength (zero
momentum transfer) limit. De-excitations of the final-state nucleus emerging
from the primary interaction are simulated using a combination of tabulated
|gamma|-ray decay schemes and an original implementation of the Hauser-Feshbach
statistical model.

Input files are provided with the code that are suitable for simulating the
charged-current process

|ve40ArCC|

coherent elastic neutrino-nucleus scattering (CEvNS) on spin-zero target
nuclei, and neutrino-electron elastic scattering on any atomic target.
Inclusion of additional reactions and targets is planned for the future.

.. |ve40ArCC| raw:: html

   <p align="center">&nu;<sub>e</sub>&nbsp;+&nbsp;<sup>40</sup>Ar&nbsp;&rarr;
   &nbsp;e<sup>&minus;</sup>&nbsp;+&nbsp;<sup>40</sup>K<sup>&ast;</sup>,</p>

.. overview-end

MARLEY follows an open-source development model and welcomes contributions of
new input files and code improvements from the community. A partial list of
potential projects for future MARLEY development is available on the developer
documentation `webpage
<http://www.marleygen.org/dev_docs.html#development-wish-list>`__.

See the `FILEMAP <FILEMAP>`__ file for a description of the full contents
of the MARLEY source code distribution.

Copyright and License
---------------------

.. copyright-start-1

.. |copy| unicode:: 0xA9 .. copyright sign

Copyright |copy| 2016-2021 Steven Gardiner gardiner@fnal.gov

MARLEY is distributed under the terms of version 3 of the `GNU General Public
License <http://www.gnu.org/licenses/gpl-3.0-standalone.html>`__ ("GPLv3") as
published by the Free Software Foundation. For the full text of that license,
please see the `COPYING <COPYING>`__ file.

.. copyright-start-2
As a matter of professional courtesy, MARLEY users are also requested to follow
the `MCnet Guidelines <https://www.montecarlonet.org/GUIDELINES>`__ for Event
Generator Authors and Users. Nevertheless, these guidelines are not legally
binding and do not limit your rights guaranteed under the GPLv3.
See the `GUIDELINES <GUIDELINES>`__ file for more details.

Citing MARLEY
-------------

.. citing-start

If you refer to MARLEY in academic work, please **always cite** the following
reference:

S. Gardiner, Simulating low-energy neutrino interactions with MARLEY,
`Comput. Phys. Commun. 269, 108123
<https://doi.org/10.1016/j.cpc.2021.108123>`__,
`arXiv:2101.11867 [nucl-th] <https://arxiv.org/abs/2101.11867>`__ (2021).

In publications which use the official reaction input files for
charged-current scattering on argon-40 (i.e., any of the files
in ``data/react`` whose names begin with ``ve40ArCC``), please also
cite the paper describing their preparation:

.. |endOfPaperTitle| raw:: html

   &nu;<sub>e</sub> scattering on <sup>40</sup>Ar,

S. Gardiner, Nuclear de-excitations in low-energy charged-current
|endOfPaperTitle|
`Phys. Rev. C 103, 044604
<https://doi.org/10.1103/PhysRevC.103.044604>`__,
`arXiv:2010.02393 [nucl-th] <https://arxiv.org/abs/2010.02393>`__ (2021).

Providing a citation for the MARLEY code itself is also encouraged and
appreciated. To maximize reproducibility of published calculations, such
citations should include the digital object identifier (DOI) associated with
the code release that was used. The DOIs for recent versions of MARLEY are
listed on the GitHub `releases webpage
<https://github.com/MARLEY-MC/marley/releases>`__ and in the right-hand column
of the Zenodo "concept DOI" `webpage
<https://doi.org/10.5281/zenodo.3901933>`__.

For convenience, recommended `BibTeX <http://www.bibtex.org/>`__ citations to
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

MARLEY is regularly `tested <https://travis-ci.org/github/MARLEY-MC/marley>`__
on both Linux and macOS platforms and is expected to work in any Unix-like
environment in which the prerequisites are installed. Building and running
MARLEY on Windows is not currently supported.

Prerequisites
~~~~~~~~~~~~~

There are three prerequisites needed to build MARLEY:

.. getting-started-end1

.. class:: open

.. getting-started-start2

.. |gte| unicode:: 0x2265 .. greater than or equal to sign

*  A C++14-compliant compiler. The following compilers are officially
   supported:

   -  `GNU Compiler Collection <https://gcc.gnu.org>`__ (GCC) |gte| 4.9.4

   -  `Clang <https://clang.llvm.org>`__ |gte| 3.5.2

*  `GNU Make <https://www.gnu.org/software/make/>`__

*  `GNU Scientific Library <https://www.gnu.org/software/gsl/>`__ (GSL)

   - MARLEY's ``Makefile`` verifies that GSL is installed by
     checking that the ``gsl-config`` script is available on the system
     ``PATH``.

On Linux machines, all three of these prerequisites will likely be available
through the standard package manager. On macOS, installing GSL may be done
using `Homebrew <https://brew.sh/>`__:

::

  brew install gsl

Although it is not required in order to build or use MARLEY, the popular `ROOT
<https://root.cern.ch>`__ data analysis framework provides convenient tools for
plotting and analyzing simulation results. Users who wish to use the optional
interface between the two codes should ensure that ROOT is installed before
building MARLEY. At build time, the optional MARLEY interface to ROOT is
enabled automatically if the ``root-config`` script is present on the system
``PATH``.

.. getting-started-end2

.. getting-started-start3

Building MARLEY
~~~~~~~~~~~~~~~

To build the code, enter the ``build/`` folder

::

    cd build

and then run GNU make

::

    make

If the build is successful, then executing

::

    ./marley

should produce the following output:

::

    Usage: marley [OPTION...] CONFIG_FILE

      -h, --help     Print this help message
      -v, --version  Print version and exit

    MARLEY home page: <http://www.marleygen.org>
    E-mail bug reports to: <support@marleygen.org>

Setting up the runtime environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``marley`` executable relies on the system environment variable ``MARLEY``
to store the full path to the root folder of the source code. This variable may
be set automatically by executing ("sourcing") the ``setup_marley.sh`` Bash
script using the ``source`` command. From within the ``build/`` folder, for
example, one may source the setup script via

::

    source ../setup_marley.sh

If generation of events is attempted without setting the ``MARLEY`` environment
variable first, then MARLEY will halt after printing the error message

::

    [ERROR]: The MARLEY environment variable is not set. Please set it (e.g.,
    by sourcing the setup_marley.sh script) and try again.

For user convenience, the ``setup_marley.sh`` script also adds the ``build/``
folder to the system ``PATH`` and to either ``LD_LIBRARY_PATH`` (Linux) or
``DYLD_LIBRARY_PATH`` (macOS).

Generating events
~~~~~~~~~~~~~~~~~

The ``marley`` executable allows the user to adjust simulation parameters via
job configuration files written in a `JSON
<https://www.json.org/json-en.html>`__-like format. The name of the
configuration file to use appears as the first (and only) command-line
argument:

::

  marley CONFIG_FILE

To generate events using an example configuration file, execute the following
command from within the ``build/`` folder after sourcing the
``setup_marley.sh`` script:

::

    marley ../examples/config/annotated.js

The program will display the MARLEY logo and diagnostic messages as it runs the
simulation. When the program terminates, a new file named ``events.ascii`` will
be present in the ``build/`` folder. This file contains the generated events
in MARLEY's native ASCII output format.

The ``annotated.js`` configuration file mentioned above is heavily commented
with explanations of the most commonly-used input parameters. Reading it serves
as a good next step for new users. When you are ready to start writing your own
configuration files, editing a copy of ``examples/config/COPY_ME.js`` is
recommended.

Full documentation for configuring MARLEY is available in section 6 of the
MARLEY `implementation paper <http://arxiv.org/abs/2101.11867>`__.

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
official webpage at http://www.marleygen.org/.

Acknowledgements
----------------

Special thanks go to

.. class:: open

- The `TALYS <http://talys.eu>`__ developers (Arjan Koning, Stéphane
  Hilaire, and Marieke Duijvestijn) for sharing their nuclear structure data

- Zero Anixter for providing an illustration of Bob Marley to be used
  in the official MARLEY logo

.. |platform| image:: https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey

.. |License: GPL v3| image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0

.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3901933.svg
   :target: https://doi.org/10.5281/zenodo.3901933

.. |Build Status| image:: https://travis-ci.org/MARLEY-MC/marley.svg?branch=main
   :target: https://travis-ci.org/MARLEY-MC/marley

.. |rel| image:: https://img.shields.io/github/v/release/MARLEY-MC/marley?include_prereleases
   :target: https://github.com/MARLEY-MC/marley/releases

.. |commits since| image:: https://img.shields.io/github/commits-since/MARLEY-MC/marley/latest/main
   :target: https://github.com/MARLEY-MC/marley/commits/main
