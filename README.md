# MARLEY (Model of Argon Reaction Low Energy Yields)

## Introduction

MARLEY (Model of Argon Reaction Low Energy Yields) is a Monte Carlo event
generator for tens-of-MeV neutrino-nucleus interactions. The current version
includes evaluated nuclear matrix elements for the charged-current reaction
<p align="center">&nu;<sub>e</sub>&nbsp;+&nbsp;<sup>40</sup>Ar&nbsp;&rarr;
&nbsp;e<sup>&minus;</sup>&nbsp;+&nbsp;<sup>40</sup>K<sup>&ast;</sup></p>
only. With the addition of new matrix elements, MARLEY could be expanded to
include more reactions and more target nuclei. Users interested in contributing
new nuclear matrix element calculations are encouraged to contact the author.

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
Detection](https://search.proquest.com/docview/2194284425), PhD thesis,
University of California, Davis (2018).

For convenience, a BibTeX citation for this thesis is given in the CITATION.bib
file.

## Getting Started

MARLEY may be built as an ordinary command-line executable, or it may be built
as product for the relocatable ups tool developed by
[Fermilab](https://fnal.gov). Those interested in doing the latter should
consult the [ups-specific
repository](https://github.com/sjgardiner/marley-ups-build). The instructions
given here describe the ordinary build procedure.

### Prerequisites
There are three prerequisites needed to build MARLEY:

- A C++14-compliant compiler. The following compilers are officially supported:
  * g++ >= 4.9.3
  * clang++ >= 3.4

- [GNU Make](https://www.gnu.org/software/make/)

- The [GNU Scientific Library](https://www.gnu.org/software/gsl/) (GSL)
  * MARLEY's Makefile verifies that GSL is installed by checking that
the `gsl-config` executable is available on the system `PATH`.

On Linux machines, all three of these prerequisites will likely be available
through the standard package manager. On Mac OSX, installing GSL may be done
using [homebrew](https://brew.sh/).

Building MARLEY on Windows has not yet been attempted. Please contact the
author (<gardiner@fnal.gov>) if you are interested in building the code
on Windows.

### Building MARLEY

To build the code, enter the build subfolder
```
cd build
```

and then run GNU make
```
make
```

If the build is successful, then executing
```
./marley --version
```

should produce the following output:
```
MARLEY (Model of Argon Reaction Low Energy Yields) 1.2.0
Copyright (C) 2016-2020 Steven Gardiner
License: GNU GPL version 3 <http://opensource.org/licenses/GPL-3.0>
This is free software: you are free to change and redistribute it.
```

### Generating events

The command-line executable allows the user to configure MARLEY via
JSON-like configuration files. To generate events using an
example configuration file, execute
```
./marley ../examples/config.js
```

The program will display the MARLEY logo and diagnostic messages as
it runs the simulation. When the program terminates, two new files will
be present in the `build/` folder:
  - `events.ascii` contains the generated events in MARLEY's native format
  - `marley.log` contains logging messages from the simulation
  
The example configuration file (`examples/config.js`) is heavily commented, and
reading it serves as a good next step for new users. When you are ready to
start writing your own configuration files, editing a copy of
`examples/COPY_ME.js` is recommended. That file contains settings similar to
those in `examples/config.js`, but it omits the long comments so that
it may be easily modified.

To use the command-line executable from outside the `build/` folder, one may optionally
run
```
sudo make install
```
to copy the `marley` executable and other files to standard system locations (e.g.,
`/usr/bin`). The installation can be removed automatically by invoking
```
sudo make uninstall
```
from the `build/` folder. For users that prefer a local installation, the `setup_marley.sh`
bash script has also been provided. Running
```
source setup_marley.sh
```
from the root MARLEY folder will add `build/` to the system `PATH`. It will also set up
a few other helpful environment variables.

## Developers

* **Steven Gardiner** - [sjgardiner](https://github.com/sjgardiner)

See also the list of
[contributors](https://github.com/sjgardiner/marley/contributors) who
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
