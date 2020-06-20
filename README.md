# MARLEY (Model of Argon Reaction Low Energy Yields)

[![Build Status](https://travis-ci.org/MARLEY-MC/marley.svg?branch=master)](https://travis-ci.org/MARLEY-MC/marley) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3901933.svg)](https://doi.org/10.5281/zenodo.3901933)

## Introduction

MARLEY (Model of Argon Reaction Low Energy Yields) is a Monte Carlo event
generator for tens-of-MeV neutrino-nucleus interactions. The current version
primarily focuses on simulations of the charged-current reaction <p
align="center">&nu;<sub>e</sub>&nbsp;+&nbsp;<sup>40</sup>Ar&nbsp;&rarr;
&nbsp;e<sup>&minus;</sup>&nbsp;+&nbsp;<sup>40</sup>K<sup>&ast;</sup>.</p>
Preparation of new reaction input files will allow MARLEY to simulate
additional reactions on more target nuclei. Users interested in extending the
functionality of MARLEY in this way are encouraged to contact the author.

## Copyright and License

Copyright &copy; 2016-2020 Steven Gardiner <gardiner@fnal.gov>

MARLEY is distributed under the terms of version 3 of the [GNU General Public
License](http://www.gnu.org/licenses/gpl-3.0-standalone.html) ("GPLv3") as
published by the Free Software Foundation. For the full text of that license,
please see the [COPYING](COPYING) file.

As a matter of professional courtesy, MARLEY users are also requested to follow
the [MCnet Guidelines](https://www.montecarlonet.org/GUIDELINES) for Event
Generator Authors and Users. Nevertheless, these guidelines are not legally
binding and do not limit your rights guaranteed under the GPLv3. See the
[GUIDELINES](GUIDELINES) file for more details.

## Citing MARLEY

If you refer to MARLEY in academic work, please __always cite__ the following
reference:

S. Gardiner, [Nuclear Effects in Neutrino
Detection](http://old.inspirehep.net/record/1802074/), PhD thesis,
University of California, Davis (2018).

Providing a citation for the MARLEY code itself is also encouraged and
appreciated. To maximize reproducibility of published calculations, such
citations should include the digital object identifier (DOI) associated with
the code release that was used. The DOIs for recent versions of MARLEY are
listed on the GitHub [releases
webpage](https://github.com/MARLEY-MC/marley/releases) and in the right-hand
column of the Zenodo "concept DOI"
[webpage](https://doi.org/10.5281/zenodo.3901933).

For convenience, recommended [BibTeX](http://www.bibtex.org/) citations to use
for the latest MARLEY release are given in the [CITATION.bib](CITATION.bib)
file.

If you use the default nuclear structure data files (strongly recommended) for
published calculations, please also give proper attribution to the developers
of the TALYS nuclear code (see
[data/structure/README.md](data/structure/README.md) for more information).

## Getting Started

MARLEY is regularly tested on both Linux and macOS platforms and is expected to
work in any Unix-like environment in which the prerequisites are installed.
Building and running MARLEY on Windows is not currently supported.

### Prerequisites
There are three prerequisites needed to build MARLEY:

- A C++14-compliant compiler. The following compilers are officially supported:

  * g++ >= 4.9.4
  * clang++ >= 3.5.2

- [GNU Make](https://www.gnu.org/software/make/)

- The [GNU Scientific Library](https://www.gnu.org/software/gsl/) (GSL)

  * MARLEY's Makefile verifies that GSL is installed by checking that
the `gsl-config` executable is available on the system `PATH`.

On Linux machines, all three of these prerequisites will likely be available
through the standard package manager. On macOS, installing GSL may be done
using [Homebrew](https://brew.sh/).

### Building MARLEY

To build the code, enter the `build/` folder
```
cd build
```
and then run GNU make
```
make
```
If the build is successful, then executing
```
./marley
```
should produce the following output:
```
Usage: marley [OPTION...] CONFIG_FILE

  -h, --help     Print this help message
  -v, --version  Print version and exit

MARLEY home page: <http://www.marleygen.org>
E-mail bug reports to: <support@marleygen.org>
```

### Configuring the runtime environment

The `marley` executable relies on the system environment variable `MARLEY`
to store the full path to the root folder of the source code. This variable
may be set automatically by executing ("sourcing") the `setup_marley.sh`
Bash script using the `source` command. From within the `build/` folder,
for example, one may source the setup script via
```
source ../setup_marley.sh
```
If generation of events is attempted without setting the `MARLEY` environment
variable first, then MARLEY will halt after printing the error message
```
[ERROR]: The MARLEY environment variable is not set. Please set it (e.g., by sourcing the setup_marley.sh script) and try again.
```
For user convenience, the `setup_marley.sh` script also adds the `build/`
folder to the system `PATH` and to either `LD_LIBRARY_PATH` (Linux) or
`DYLD_LIBRARY_PATH` (macOS).

### Generating events

The `marley` command-line executable allows the user to adjust simulation
parameters via JSON-like configuration files. To generate events using an
example configuration file, execute the following command from within the
`build/` folder after sourcing the `setup_marley.sh` script:
```
marley ../examples/config/annotated.js
```

The program will display the MARLEY logo and diagnostic messages as
it runs the simulation. When the program terminates, two new files will
be present in the `build/` folder:
  - `events.ascii` contains the generated events in MARLEY's native ASCII format
  - `marley.log` contains logging messages from the simulation
  
The example configuration file (`examples/config.js`) is heavily commented, and
reading it serves as a good next step for new users. When you are ready to
start writing your own configuration files, editing a copy of
`examples/COPY_ME.js` is recommended.

## Developers

* **Steven Gardiner** - [sjgardiner](https://github.com/sjgardiner)

See also the list of
[contributors](https://github.com/MARLEY-MC/marley/contributors) who
participated in this project.

## Website

[Doxygen](https://www.doxygen.org) documentation for the latest version of
MARLEY may be found on the official webpage at <http://www.marleygen.org/>.

## Acknowledgements

Special thanks go to

* The [TALYS](http://talys.eu) developers (Arjan Koning, Stéphane Hilaire, and
  Marieke Duijvestijn) for sharing their nuclear structure data

* Zero Anixter for providing an illustration of Bob Marley to be used in the
  official MARLEY logo
