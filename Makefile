# Top-level Makefile for MARLEY
# Dispatches to CMake (preferred) or the hand-written GNU Make recipe (make/build.mk).
#
# Usage:
#   make                   — auto-detect cmake; use it if found
#   make IGNORE_CMAKE=1    — force use of the hand-written make/build.mk recipe
#   make IGNORE_ROOT=1     — pass through to either build system
#   make IGNORE_HEPMC3=1   — pass through to either build system
#   make IGNORE_GSL=1      — use built-in GSL subset (pass through to either build system)
#   make debug             — debug build
#   make clean             — remove build/ entirely
#   make install           — install under $(prefix) (default /usr/local)
#   make help              — show this message
#
# Install locations are given with the GNU names on both build paths:
#   make install prefix=/opt/marley libdir=lib64 DESTDIR=/tmp/stage

CMAKE := $(shell command -v cmake 2>/dev/null)

USE_CMAKE :=
ifndef IGNORE_CMAKE
ifdef CMAKE
USE_CMAKE := 1
endif
endif

CMAKE_FLAGS :=
ifdef IGNORE_ROOT
CMAKE_FLAGS += -DIGNORE_ROOT=ON
endif
ifdef IGNORE_HEPMC3
CMAKE_FLAGS += -DIGNORE_HEPMC3=ON
endif
ifdef IGNORE_GSL
CMAKE_FLAGS += -DIGNORE_GSL=ON
endif

# The install-location variables are spelled with the GNU names on the command
# line, because that is what packaging recipes already pass. CMake spells them
# differently, so this Makefile -- the single front door over both build paths
# -- translates rather than leaving 'make install prefix=...' to mean two
# different things depending on whether cmake happens to be installed.
#
# These are configure-time settings for CMake, not install-time ones: libdir in
# particular also determines where the build tree puts libMARLEY (the build
# tree mirrors the install layout), and the relative form is baked into the
# RUNPATH. A cached value must therefore not win over one given here, so any of
# them forces a reconfigure below.
CMAKE_DIR_FLAGS :=
ifdef prefix
CMAKE_DIR_FLAGS += -DCMAKE_INSTALL_PREFIX='$(prefix)'
endif
ifdef libdir
CMAKE_DIR_FLAGS += -DCMAKE_INSTALL_LIBDIR='$(libdir)'
endif
ifdef bindir
CMAKE_DIR_FLAGS += -DCMAKE_INSTALL_BINDIR='$(bindir)'
endif
ifdef includedir
CMAKE_DIR_FLAGS += -DCMAKE_INSTALL_INCLUDEDIR='$(includedir)'
endif

# Reconfigure only when it would change something: an unconditional 'cmake -S'
# on every build is slow, but a stale cache silently overriding an explicit
# command-line prefix is the defect being repaired here.
ifeq ($(strip $(CMAKE_DIR_FLAGS)),)
CMAKE_CONFIGURE = @test -d build/CMakeFiles || cmake -S . -B build $(CMAKE_FLAGS)
else
CMAKE_CONFIGURE = cmake -S . -B build $(CMAKE_FLAGS) $(CMAKE_DIR_FLAGS)
endif

# 'cmake --install' has no --destdir flag; the environment variable is the
# documented interface, and it is read at install time only.
CMAKE_INSTALL_ENV :=
ifdef DESTDIR
CMAKE_INSTALL_ENV := DESTDIR='$(DESTDIR)'
endif

PASSTHROUGH_VARS :=
ifdef IGNORE_ROOT
PASSTHROUGH_VARS += IGNORE_ROOT=$(IGNORE_ROOT)
endif
ifdef IGNORE_HEPMC3
PASSTHROUGH_VARS += IGNORE_HEPMC3=$(IGNORE_HEPMC3)
endif
ifdef IGNORE_GSL
PASSTHROUGH_VARS += IGNORE_GSL=$(IGNORE_GSL)
endif
ifdef CXX
PASSTHROUGH_VARS += CXX='$(CXX)'
endif
ifdef CXXFLAGS
PASSTHROUGH_VARS += CXXFLAGS='$(CXXFLAGS)'
endif
ifdef prefix
PASSTHROUGH_VARS += prefix='$(prefix)'
endif
ifdef libdir
PASSTHROUGH_VARS += libdir='$(libdir)'
endif
ifdef bindir
PASSTHROUGH_VARS += bindir='$(bindir)'
endif
ifdef includedir
PASSTHROUGH_VARS += includedir='$(includedir)'
endif
ifdef DESTDIR
PASSTHROUGH_VARS += DESTDIR='$(DESTDIR)'
endif

.PHONY: all clean debug reconfigure test docs doxygen install uninstall help marg4

all:
ifdef USE_CMAKE
	$(CMAKE_CONFIGURE)
	$(MAKE) -C build
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) $(PASSTHROUGH_VARS)
endif

debug:
ifdef USE_CMAKE
	cmake -S . -B build $(CMAKE_FLAGS) $(CMAKE_DIR_FLAGS) -DCMAKE_BUILD_TYPE=Debug
	$(MAKE) -C build
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) debug $(PASSTHROUGH_VARS)
endif

test:
ifdef USE_CMAKE
	$(CMAKE_CONFIGURE)
	$(MAKE) -C build test
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) test $(PASSTHROUGH_VARS)
endif

marg4:
ifdef USE_CMAKE
	$(CMAKE_CONFIGURE)
	cmake --build build --target marg4
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) marg4 $(PASSTHROUGH_VARS)
endif

docs:
ifdef USE_CMAKE
	$(CMAKE_CONFIGURE)
	cmake --build build --target docs $(CMAKE_BUILD_FLAGS)
	$(MAKE) -C build docs
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) docs $(PASSTHROUGH_VARS)
endif

doxygen:
ifdef USE_CMAKE
	$(CMAKE_CONFIGURE)
	cmake --build build --target doxygen $(CMAKE_BUILD_FLAGS)
	$(MAKE) -C build doxygen
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) doxygen $(PASSTHROUGH_VARS)
endif

install:
ifdef USE_CMAKE
	$(CMAKE_CONFIGURE)
# Build before installing: an install-location variable given here may have
# just forced a reconfigure, and libdir in particular moves the build tree's
# own library directory, so the artifacts to install may not exist yet.
	$(MAKE) -C build
	$(CMAKE_INSTALL_ENV) cmake --install build
else
	@mkdir -p build
	$(MAKE) -C build -f $(CURDIR)/make/build.mk TOP_DIR=$(CURDIR) install $(PASSTHROUGH_VARS)
endif

# Both build paths dispatch to the same script, which reads the receipt the
# installer left in the prefix. It needs no build tree and no CMake cache --
# only a prefix -- so any checkout can uninstall any installation.
uninstall:
	src/scripts/uninstall-marley --prefix "$(DESTDIR)$(or $(prefix),/usr/local)"

reconfigure:
ifdef USE_CMAKE
	$(RM) -rf build/CMakeFiles
	cmake -S . -B build $(CMAKE_FLAGS) $(CMAKE_DIR_FLAGS)
else
	@echo "reconfigure is only supported for the CMake build path."
endif

clean:
	$(RM) -rf build

help:
	@echo "Usage:"
	@echo "  make                   Auto-detect cmake and use it if found"
	@echo "  make IGNORE_CMAKE=1    Force GNU Make path (make/build.mk)"
	@echo "  make IGNORE_ROOT=1     Ignore ROOT"
	@echo "  make IGNORE_HEPMC3=1   Ignore external HepMC3"
	@echo "  make IGNORE_GSL=1      Use built-in GSL subset"
	@echo "  make marg4             Build the MARLEY+Geant4 example (if Geant4 is available)"
	@echo "  make [debug|test|docs|doxygen|install|uninstall|clean|reconfigure]"
	@echo ""
	@echo "  'make uninstall' asks for confirmation, and refuses rather than"
	@echo "  guessing if it cannot. For a script or CI, call the uninstaller"
	@echo "  directly: src/scripts/uninstall-marley --prefix DIR --yes"
	@echo ""
	@echo "Install locations (GNU names, honoured on both build paths):"
	@echo "  prefix=DIR             Installation prefix (default /usr/local)"
	@echo "  libdir=NAME            Library directory name (default lib)"
	@echo "  bindir=NAME            Executable directory name (default bin)"
	@echo "  includedir=NAME        Header directory name (default include)"
	@echo "  DESTDIR=DIR            Staging directory prepended to every path"
