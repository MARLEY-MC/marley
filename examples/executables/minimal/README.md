# "Minimal" executable examples

For a variety of applications, it may be helpful to create external C++
programs that make use of one or more MARLEY classes. To aid users in building
such programs, a Bash script called *marley-config* is generated and placed in
the build/ folder whenever MARLEY is successfully built. The three programs
in this folder may each be compiled using the marley-config script by executing
a command of the form

```
g++ -o EXECUTABLE_NAME $(marley-config --cflags --libs) SOURCE_NAME
```

where EXECUTABLE_NAME is the file name for the compiled executable and
SOURCE_NAME is the name of a C++ source file in this folder. The syntax shown
above assumes that the build/ folder has previously been added to the system
PATH (e.g., by sourcing the setup_marley.sh script in the root source code
folder).

The three example programs provided in this folder are

* efr.cc: Reads a single MARLEY output file specified as the first command-line
  argument. The flux-averaged total cross section is printed to stdout,
  followed by each event in the ASCII output format.

* evgen.cc: Generates 10 MARLEY events using the job configuration file
  /home/config.js. This file name is currently hard-coded. Each event is
  printed to stdout in ASCII format after being generated.

* mass_40Ar.cc: Prints a brief message to stdout with two pieces of
  information: (1) whether or not MARLEY was built with ROOT support
  enabled, and (2) the atomic mass of 40Ar in MeV/c^2

Examples of using the marley-config script in a Makefile may be seen in
examples/executables/build/Makefile and examples/marg4/build/Makefile.
