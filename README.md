# MARLEY (Model of Argon Reaction Low Energy Yields)

## Introduction

MARLEY (Model of Argon Reaction Low Energy Yields) is a Monte Carlo event
generator for tens-of-MeV neutrino-nucleus interactions in liquid argon. The
current version includes matrix element evaluations for the charged-current
reaction
&nu;<sub>e</sub>&nbsp;+&nbsp;<sup>40</sup>Ar&nbsp;&rarr;&nbsp;e<sup>&minus;
</sup>&nbsp;+&nbsp;<sup>40</sup>K<sup>*</sup> only.
With some additional matrix element evaluations, MARLEY could be expanded to
include more reactions and more target nuclei. Users interested in creating
reaction matrix element evaluations for new reactions or new nuclear targets
are encouraged to contact the MARLEY development team
(<support@marleygen.org>).

## Licensing

MARLEY is distributed under the terms of the [BSD 2-clause
license](http://opensource.org/licenses/BSD-2-Clause). Please see the LICENSE
file for details.

## Getting Started

MARLEY may be built as an ordinary command-line executable, or it may be built
as product for the relocatable ups tool developed by
[Fermilab](https://fnal.gov). Those interested in doing the latter should
consult the [ups-specific
repository](https://github.com/sjgardiner/marley-ups-build). The instructions
given here describe the ordinary build procedure.

### Prerequisites
There are only two prerequisites needed to build MARLEY:

- A C++14-compliant compiler. The following compilers are officially supported:
  * g++ >= 4.9
  * clang++ >= 3.4

- [GNU Make](https://www.gnu.org/software/make/)

Building MARLEY on Windows has not yet been attempted. Please contact the
developers (<support@marleygen.org>) if you are interested in building the code
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
MARLEY (Model of Argon Reaction Low Energy Yields) 1.1.0
Copyright (C) 2016-2018 Steven Gardiner
License BSD2C: BSD 2-Clause <http://opensource.org/licenses/BSD-2-Clause>
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

## Authors

* **Steven Gardiner** - [sjgardiner](https://github.com/sjgardiner)

Non-developer MARLEY authors:
**Christopher Grant**, **Emilija Pantic**, **Robert Svoboda**

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
