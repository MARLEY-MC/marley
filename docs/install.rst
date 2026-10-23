============
Installation
============

All MARLEY features are available while running from within the build tree, so
installation is optional. It is recommended for users who want to package
MARLEY with other software or run the code from a more standard location.

Syntax
------

To install MARLEY under the default prefix ``/usr/local``, run

::

  make install

from the root folder of the source code. The ``prefix`` variable may be used
to specify a different install location:

::

  make install prefix=/opt/marley

The following variables are accepted by the top-level Makefile and applied
consistently for both CMake and GNU Make builds:

``prefix``
  Root of the installation. Defaults to ``/usr/local``.

``libdir``
  Name of the library directory. Defaults to ``lib``.

``bindir``
  Name of the executable directory. Defaults to ``bin``.

``includedir``
  Name of the header file directory. Defaults to ``include``.

``DESTDIR``
  Staging directory prepended to every installed path. Used by packaging
  recipes, which install into a temporary tree before building an archive.
  It does not affect the paths recorded inside the installed files.

Two restrictions apply to the directory names above. Each must refer to an
immediate sub-folder of the prefix, so ``libdir=lib64`` is accepted but
``libdir=lib/x86_64-linux-gnu`` is not. This restriction is imposed to ensure
that MARLEY can automatically find the installed data files. The data directory
name is fixed at ``share`` for the same reason.

``libdir`` is also a build-time setting rather than an install-time one,
because the relative path it names is recorded in the executables when they
are linked. If specified during the build, it should also be given with the
same value during installation:

::

  make libdir=lib64
  make install libdir=lib64

Layout of the installed files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

  <prefix>/bin/marley                      the command-line executable
  <prefix>/bin/marley-config               reports MARLEY build flags
  <prefix>/bin/setup_marley.sh             optional environment script
  <prefix>/bin/mroot                       ROOT with MARLEY loaded (ROOT builds)
  <prefix>/<libdir>/libMARLEY.so           the MARLEY shared library
  <prefix>/<libdir>/pkgconfig/marley.pc    pkg-config metadata
  <prefix>/<libdir>/cmake/MARLEY/          CMake package (CMake builds only)
  <prefix>/<includedir>/marley/            the public MARLEY header files
  <prefix>/share/marley/data/              input physics data and configuration
  <prefix>/share/marley/examples/          example configurations and programs
  <prefix>/share/marley/install-receipt    manifest used by the uninstaller

Running after installation
--------------------------

MARLEY locates its own data files relative to the shared library's position on
disk, and the executables record where to find the library, so running directly
from an installed folder, e.g.,

::

  /opt/marley/bin/marley my_config.js

is expected to work without setting any environment variables. To avoid typing
the path to the executable, one may optionally add ``<prefix>/bin`` to ``PATH``.

The script ``<prefix>/bin/setup_marley.sh`` is installed to support (unusual)
cases that require changes to the system environment. An example is using
MARLEY classes from within a Jupyter notebook running the ROOT C++ kernel.
Sourcing this script will modify several environment variables to account for
MARLEY's install location:

::

  source /opt/marley/bin/setup_marley.sh

A MARLEY installation is relocatable. That is, the installation prefix may be
moved or copied elsewhere and is expected to continue to work:

::

  mv /opt/marley /opt/marley-2.1.0

To report the directories MARLEY uses to find data files, use the command

::

  marley info

Additional directories may be prepended to the search order through the
``MARLEY_SEARCH_PATH`` environment variable, which holds a colon-separated list.
Directories named there are searched before the installed data area, allowing
user-supplied input files to take precedence over the official ones.

Building other programs against MARLEY
--------------------------------------

Three tools are provided for external codes that interface with MARLEY.

**marley-config** prints the flags used to build MARLEY. This script locates
the installation from its own position, so no arguments are needed:

::

  g++ -o my_program my_program.cc $(marley-config --cflags --libs)

Valid options also include ``--prefix``, ``--bindir``, ``--libdir``,
``--includedir``, ``--cxx-std``, ``--use-root``, ``--version`` and
``--git-revision``. An example Makefile using ``marley-config`` may
be found at ``<prefix>/share/marley/examples/marg4/Makefile``.

**pkg-config** reads ``marley.pc`` and behaves similarly. It may be preferred
for some projects. Point ``PKG_CONFIG_PATH`` at the installation if it is not
in a default location:

::

  export PKG_CONFIG_PATH=/opt/marley/lib/pkgconfig:${PKG_CONFIG_PATH}
  g++ -o my_program my_program.cc $(pkg-config --cflags --libs marley)

**CMake** projects should typically use the exported package, which is
installed only when building MARLEY with CMake:

::

  find_package(MARLEY 2.1 REQUIRED)
  target_link_libraries(my_program PRIVATE MARLEY::MARLEY)

Add the installation prefix to ``CMAKE_PREFIX_PATH`` if it is outside the
default search path. The imported target carries the include directories, the
C++ standard used to compile MARLEY, and its public dependencies.

Uninstalling
------------

The installation produces a manifest ``<prefix>/share/marley/install-receipt``
that lists everything that was placed in the prefix. From the root folder
of the MARLEY source code,

::

  make uninstall prefix=/opt/marley

reads that manifest and removes the installation, asking for confirmation
first. One may call the uninstaller directly to skip the confirmation prompt:

::

  src/scripts/uninstall-marley --prefix /opt/marley --yes

The uninstaller is not itself installed. Since the manifest lists the content to
be removed, any MARLEY source tree may be used to perform the removal. The files
listed in the manifest may also be deleted by hand.
