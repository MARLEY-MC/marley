# Names of various MARLEY directories
TOP_DIR ?= $(CURDIR)
BUILD_DIR = $(TOP_DIR)/build
INCLUDE_DIR = $(TOP_DIR)/include
SRC_DIR = $(TOP_DIR)/src

# ---------------------------------------------------------------------------
# Install directories.
#
# These are defined here, rather than next to the install target, because the
# build tree mirrors the install layout: the library is built into
# $(BUILD_DIR)/$(libdir_rel), and SHARED_LIB below needs that value.
#
# libdir, bindir and includedir must resolve to exactly ONE directory below
# the prefix, because the data root is found by a single ".." ascent from the
# library's own location (design spec 2.2, 4.3). Both the GNU-conventional
# absolute spelling and a bare segment are accepted:
#
#   make install libdir=lib64              ->  lib64     OK
#   make install libdir=/usr/lib64         ->  lib64     OK  (prefix=/usr)
#   make install libdir=/usr/lib/x86_64-…  ->  ERROR     multi-segment
#   make install libdir=/opt/other/lib     ->  ERROR     not under the prefix
#
# NOTE: libdir is a BUILD-TIME setting, not an install-time one. It is baked
# into the RUNPATH of the executables at link time, so it must be given to
# both invocations: make libdir=lib64 && make install libdir=lib64. (This was
# already true before it was written down -- the old recipe put $(libdir) on
# the link line -- which is why "make && make install libdir=lib64" has never
# produced a working install.)
# ---------------------------------------------------------------------------

# Strip any trailing slash so that the patsubst below matches. "override" is
# required throughout this block: a variable set on the command line beats an
# ordinary makefile assignment, so without it every normalisation below would
# be silently skipped for precisely the invocations it exists to handle.
override prefix := $(patsubst %/,%,$(or $(prefix),/usr/local))
override exec_prefix := $(patsubst %/,%,$(or $(exec_prefix),$(prefix)))

bindir ?= $(exec_prefix)/bin
libdir ?= $(exec_prefix)/lib
includedir ?= $(prefix)/include
# datadir is fixed: unlike the others it is used for DESCENT, so its name
# must be known at runtime rather than merely counted (design spec 4.3).
datadir = $(prefix)/share

bindir_rel := $(patsubst $(prefix)/%,%,$(bindir))
libdir_rel := $(patsubst $(prefix)/%,%,$(libdir))
includedir_rel := $(patsubst $(prefix)/%,%,$(includedir))

define check_install_dir
ifeq ($$(strip $$($(1)_rel)),)
$$(error $(1) must name one directory below prefix=$$(prefix))
endif
ifneq (,$$(findstring /,$$($(1)_rel)))
$$(error $(1)=$$($(1)) must be exactly one directory below prefix=$$(prefix) (got $$($(1)_rel)); use e.g. $(1)=lib64)
endif
endef

$(eval $(call check_install_dir,bindir))
$(eval $(call check_install_dir,libdir))
$(eval $(call check_install_dir,includedir))

# Re-normalise so that a bare segment ("libdir=lib64") still yields an
# absolute install destination. Without this, $(DESTDIR)$(libdir) would be
# a path relative to the build directory.
override bindir := $(prefix)/$(bindir_rel)
override libdir := $(prefix)/$(libdir_rel)
override includedir := $(prefix)/$(includedir_rel)

# The build tree mirrors the install layout (design spec 2.2)
BUILD_LIBDIR = $(BUILD_DIR)/$(libdir_rel)

# Default version of the C++ standard to use for compilation. The effective
# value is computed once ROOT has been found; see "the effective C++
# standard" below. Whether the user asked for one explicitly has to be
# recorded here, while the variable's origin is still intact -- the later
# "override" changes it.
CXX_STD=c++17
ifeq ($(origin CXX_STD),command line)
  CXX_STD_EXPLICIT := yes
else ifeq ($(origin CXX_STD),environment)
  CXX_STD_EXPLICIT := yes
else
  CXX_STD_EXPLICIT :=
endif

UNKNOWN_REV=unknown version

# If the target is "debug", then include debugging information and build the
# marley executable with optimization turned off. Prepend the default -std
# option so that it will be overridden if the user has manually specified one
# in CXXFLAGS.
ifeq ($(MAKECMDGOALS),debug)
  MARLEY_LOG_LEVEL := TRACE
  override CXXFLAGS := -std=$(CXX_STD) $(CXXFLAGS) -O0 -g \
    -DMARLEY_COMPILED_LOG_LEVEL=$(MARLEY_LOG_LEVEL)
else
  # Otherwise, use full optimization and do not include debugging info. Prepend
  # the default -O3 and -std options here in case the user wants to manually
  # override them in CXXFLAGS.
  MARLEY_LOG_LEVEL := INFO
  override CXXFLAGS := -O3 -std=$(CXX_STD) $(CXXFLAGS) \
    -DMARLEY_COMPILED_LOG_LEVEL=$(MARLEY_LOG_LEVEL)
endif

UNAME_S := $(shell uname -s)
# DL_LIBS: dladdr() is used for self-location (design spec 2.1). On glibc
# older than 2.34 it lives in libdl and needs the flag; from 2.34 libdl was
# folded into libc and the flag is harmless. MARLEY targets Scientific
# Linux 7 (glibc 2.17), so it is passed unconditionally on Linux. On macOS
# dladdr is in libSystem, which is always linked.
#
# RPATH_LIB / RPATH_EXE: every linked artifact locates its dependencies
# relative to its own position, so no absolute path is baked into anything
# and the installed tree can be moved (design spec 4.2).
#
# libMARLEY needs its own RUNPATH, and this is a correctness requirement
# rather than a convenience: ld.so(8) states that DT_RUNPATH directories
# "do not apply to those objects' children", so the executable's RUNPATH
# cannot help libMARLEY find libMARLEYHepMC3.
#
# --enable-new-dtags asks for DT_RUNPATH rather than DT_RPATH. Which one
# -Wl,-rpath produces is a binutils configure-time default, and the whole
# reason for preferring RUNPATH -- that a user's LD_LIBRARY_PATH still wins
# -- depends on getting it, so it is requested explicitly.
#
# The '$$ORIGIN' quoting is deliberate: $$ escapes the dollar from Make, and
# the single quotes keep it from the shell, so the linker sees $ORIGIN.
ifeq ($(UNAME_S),Darwin)
  SHARED_LIB_SUFFIX=dylib
  DL_LIBS :=
  NEW_DTAGS :=
  RPATH_LIB := -Wl,-rpath,@loader_path
  RPATH_EXE = -Wl,-rpath,@loader_path/../$(libdir_rel)
  # Without this the recorded install name is the build-tree path, and the
  # library is not relocatable at all
  INSTALL_NAME_FLAG = -Wl,-install_name,@rpath/$(SHARED_LIB_FILE)
  HEPMC3_INSTALL_NAME_FLAG = -Wl,-install_name,@rpath/libMARLEYHepMC3.$(SHARED_LIB_SUFFIX)
else ifeq ($(UNAME_S),Linux)
  SHARED_LIB_SUFFIX=so
  DL_LIBS := -ldl
  NEW_DTAGS := -Wl,--enable-new-dtags
  RPATH_LIB := -Wl,-rpath,'$$ORIGIN'
  RPATH_EXE = -Wl,-rpath,'$$ORIGIN/../$(libdir_rel)'
  INSTALL_NAME_FLAG =
  HEPMC3_INSTALL_NAME_FLAG =
else
  $(warning Unrecognized operating system encountered.)
  SHARED_LIB_SUFFIX=so
  DL_LIBS := -ldl
  NEW_DTAGS := -Wl,--enable-new-dtags
  RPATH_LIB := -Wl,-rpath,'$$ORIGIN'
  RPATH_EXE = -Wl,-rpath,'$$ORIGIN/../$(libdir_rel)'
  INSTALL_NAME_FLAG =
  HEPMC3_INSTALL_NAME_FLAG =
endif

# EXT_RPATH: the same thing for dependencies MARLEY does NOT ship. $ORIGIN
# above locates MARLEY's own libraries; it says nothing about where ROOT, an
# external HepMC3 or GSL live, and those are not relocatable with us.
#
# Without this an installed MARLEY built against a ROOT outside the loader's
# default directories cannot start at all -- it links fine, installs fine,
# and then dies before main() with "libCore.so: cannot open shared object
# file". macOS hides the whole problem, because ROOT's dylibs record absolute
# install names and need no RUNPATH; this is a Linux-only failure, and it is
# why it survived local testing. The CMake build has the same gap for the
# same reason, closed there by CMAKE_INSTALL_RPATH_USE_LINK_PATH plus an
# equivalent explicit harvest.
#
# A fallback, not an override: NEW_DTAGS requests DT_RUNPATH, which ld.so
# searches after LD_LIBRARY_PATH, so a user's own thisroot.sh still wins.
#
# Recursively expanded (=, not :=) because ROOT_LDFLAGS and friends are set
# later, inside conditionals. Only absolute directories are taken: a relative
# -L would become a RUNPATH resolved against the current directory at load
# time, which is both useless and a hazard. Directories the loader already
# searches are dropped so that installed binaries do not look machine-
# specific for no gain.
comma := ,
SYSTEM_LIBDIRS := /lib /lib64 /usr/lib /usr/lib64 $(wildcard /usr/lib/*-linux-gnu)
EXT_LDFLAGS_ALL = $(GSL_LDFLAGS) $(ROOT_LDFLAGS) $(HEPMC3_LDFLAGS)
EXT_LIBDIRS_RAW = \
  $(patsubst -L%,%,$(filter -L%,$(EXT_LDFLAGS_ALL))) \
  $(patsubst -Wl$(comma)-rpath$(comma)%,%,\
    $(filter -Wl$(comma)-rpath$(comma)%,$(EXT_LDFLAGS_ALL)))
EXT_LIBDIRS = $(filter-out $(SYSTEM_LIBDIRS) $(abspath $(BUILD_LIBDIR)),\
  $(sort $(filter /%,$(EXT_LIBDIRS_RAW))))
EXT_RPATH = $(foreach d,$(EXT_LIBDIRS),-Wl,-rpath,$(d))

# ldconfig refreshes the Linux shared-library cache. It does not exist on
# macOS, where calling it unguarded made 'make install' fail outright, and it
# is meaningless when staging into DESTDIR, since nothing is at its final
# location yet. It is also a no-op for any prefix the cache does not cover and
# an error for an unprivileged user, so it is never fatal: refreshing a cache
# is not part of installing the files.
ifeq ($(UNAME_S),Darwin)
  LDCONFIG_CMD := true
else
  LDCONFIG_CMD := $(if $(DESTDIR),true,ldconfig 2> /dev/null || true)
endif

# If available, use git to determine the hash of the current version
# of MARLEY. If not, set it to the string "unknown version"

# Check if git is available on the system path.
GIT := $(shell command -v git 2> /dev/null)

ifndef GIT
  GIT_REVISION := $(UNKNOWN_REV)
else
  # Also check that we're inside of a folder managed by a git repository. If
  # we're not, then git rev-parse won't work
  GIT_REVPARSE_CODE := $(shell git rev-parse 2> /dev/null && echo "$$?")
  ifeq ($(GIT_REVPARSE_CODE),0)

    GIT_REVISION := $(shell git rev-parse --short HEAD)

    # Verify that the working tree and index are identical to HEAD. If they
    # are not (because there are uncommitted changes to at least one tracked
    # file), then append "-dirty" to the git revision variable.
    GIT_DIFF_INDEX_CODE := $(shell git diff-index --quiet HEAD \
      2> /dev/null && echo "$$?")

    ifneq ($(GIT_DIFF_INDEX_CODE),0)
      GIT_REVISION := $(GIT_REVISION)-dirty
    endif

  else
    GIT_REVISION := $(UNKNOWN_REV)
  endif

endif

# If the .VERSION file exists, use its contents as the
# MARLEY version number. The .VERSION file should only
# be added to the source tree in tagged releases.
ifneq (,$(wildcard $(TOP_DIR)/.VERSION))
  MARLEY_VERSION := $(shell cat $(TOP_DIR)/.VERSION)

  VERSION_PREFIX="v"

else

  # Also use the git revision as the version number in this case
  MARLEY_VERSION := $(GIT_REVISION)

endif

# Define the MARLEY_VERSION preprocessor macro
override CXXFLAGS += -DMARLEY_VERSION="\"$(MARLEY_VERSION)\""

# Define the MARLEY_GIT_REVISION preprocessor macro
override CXXFLAGS += -DMARLEY_GIT_REVISION="\"$(GIT_REVISION)\""

SHARED_LIB_NAME := MARLEY
SHARED_LIB_FILE := lib$(SHARED_LIB_NAME).$(SHARED_LIB_SUFFIX)
SHARED_LIB := $(BUILD_LIBDIR)/$(SHARED_LIB_FILE)

TEST_EXECUTABLE = $(BUILD_DIR)/bin/martest
TEST_OBJECTS = $(notdir $(patsubst %.cc,%.o,$(wildcard $(SRC_DIR)/tests/*.cc)))

all: marley marley-pc setup-script
debug: marley
test: $(TEST_EXECUTABLE)

# Skip lots of initialization if all we want is "make clean/uninstall"
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),uninstall)

  # Use g++ as the default compiler
  CXX ?= g++
  override CXXFLAGS += -I$(INCLUDE_DIR) -Wall -Wextra -Wpedantic
  override CXXFLAGS += -Wno-error=unused-parameter -Wcast-align

  # Add extra compiler flags for recognized compilers (currently just gcc
  # and clang)
  CXXVERSION = $(shell $(CXX) --version)
  COMPILER_VERSION := $(word 3, $(CXXVERSION))
  ifneq (,$(findstring clang,$(CXXVERSION)))
    # clang
    $(info Compiling using version $(COMPILER_VERSION) of clang)

    # The ROOT headers trigger clang's no-keyword-macro warning, so
    # disable it. Also disable (for now) warnings for braces around
    # initialization of subobjects (overkill in the meta_numerics header)
    CXXFLAGS += -Wno-keyword-macro -Wno-missing-braces
  else
    ifneq (,$(or $(findstring GCC,$(CXXVERSION)), $(findstring g++,$(CXXVERSION))))
      # gcc
      $(info Compiling using version $(COMPILER_VERSION) of GCC)
      # Linking to ROOT libraries can be problematic on distributions (e.g.,
      # Ubuntu) that set the g++ flag -Wl,--as-needed by default (see
      # http://www.bnikolic.co.uk/blog/gnu-ld-as-needed.html for details), so
      # disable this behavior on Linux.
      ifneq ($(UNAME_S),Darwin)
        override CXXFLAGS += -Wl,--no-as-needed
      endif
    endif
  endif

  OBJECTS := $(notdir $(patsubst %.cc,%.o,$(wildcard $(SRC_DIR)/*.cc $(SRC_DIR)/app/*.cc)))
  OBJECTS := $(filter-out marley.o marley_root.o, $(OBJECTS))
  OBJECTS := $(filter-out OutputFileRoot.o, $(OBJECTS))
  OBJECTS := $(filter-out marley_hepmc3.o, $(OBJECTS))

  # Get information about the GNU Scientific Library installation
  ifndef IGNORE_GSL
    GSLCONFIG := $(shell which gsl-config)
    ifeq (, $(GSLCONFIG))
      $(info Could not find a working GNU Scientific Library installation.)
      $(info A built-in GSL library will be used instead.)
    else
      # Add the appropriate compilation flags to use the system GNU Scientific
      # Library
      GSL_CXXFLAGS := $(shell $(GSLCONFIG) --cflags) -DMARLEY_FOUND_GSL
      GSL_LDFLAGS := $(shell $(GSLCONFIG) --libs)
      FOUND_GSL := TRUE
    endif
  else
    $(info Ignoring any GSL installations that may be present.)
  endif

  ifndef FOUND_GSL
    # Use the built-in GSL subset
    GSL_CXXFLAGS :=
    GSL_LDFLAGS := -lm

    # Keep built-in GSL symbols hidden to avoid clashes with any system GSL
    # loaded by other libraries in the host environment.
    marley_gsl.o: CXXFLAGS += -fvisibility=hidden

    # The MARLEY shared library doesn't expose the built-in GSL symbols, so
    # we will need to link them into the test executable directly. In order
    # to expose the GSL functions, we need to compile this version of the
    # built-in GSL object file *without* -fvisibility=hidden. The special
    # sufixx .test.o uses a rule below to recompile the same source into
    # a different object
    TEST_OBJECTS += marley_gsl.test.o

  endif

  # Get information about the HepMC3 installation
  ifndef IGNORE_HEPMC3
    HEPMC3CONFIG := $(shell which HepMC3-config)
    ifeq (, $(HEPMC3CONFIG))
      $(info Could not find a working HepMC3 installation.)
      $(info A built-in HepMC3 library will be used instead.)
    else
      # Add the appropriate compilation flags to use the system HepMC3
      # implementation
      HEPMC3_INCDIR   := $(shell $(HEPMC3CONFIG) --includedir)
      HEPMC3_CXXFLAGS := $(shell $(HEPMC3CONFIG) --cflags) -DMARLEY_FOUND_HEPMC3
      HEPMC3_LDFLAGS  := $(shell $(HEPMC3CONFIG) --libs)
      # Take nothing from HepMC3-config on trust. A broken script (the
      # Homebrew one fails without GNU coreutils' greadlink, for instance)
      # succeeds, prints nothing useful, and leaves an empty include
      # directory -- which surfaces thirty lines later as
      # "'HepMC3/FourVector.h' file not found", naming neither the script
      # nor the empty variable. CMake's find_package fails loudly here;
      # this makes the Make path do the same.
      ifeq (,$(wildcard $(HEPMC3_INCDIR)/HepMC3))
        $(error HepMC3-config ($(HEPMC3CONFIG)) did not report a usable \
include directory (got "$(HEPMC3_INCDIR)"). Check that it runs correctly, \
or build against the bundled HepMC3 subset with IGNORE_HEPMC3=1)
      endif
      FOUND_HEPMC3    := TRUE
    endif
  else
    $(info Ignoring any HepMC3 installations that may be present.)
  endif

  # The flags MARLEY publishes for an external HepMC3. MARLEY_FOUND_HEPMC3 is
  # an internal build switch -- no installed header tests it -- so it must not
  # reach a consumer's compiler (design spec 6.2). Empty on a vendored build,
  # where the include root is MARLEY's own and the body derives it.
  HEPMC3_PUBLIC_CFLAGS = $(if $(FOUND_HEPMC3),\
    $(filter-out -DMARLEY_FOUND_HEPMC3,$(HEPMC3_CXXFLAGS)))
  HEPMC3_PUBLIC_LIBS = $(if $(FOUND_HEPMC3),$(HEPMC3_LDFLAGS))

  ifndef FOUND_HEPMC3
    # Use the built-in HepMC3 subset (implementation in src/marley_hepmc3.cc,
    # headers in include/builtin/HepMC3/).
    #
    # IMPORTANT: The -I flag below must remain conditional on FOUND_HEPMC3
    # being unset. Making it unconditional would expose include/builtin/HepMC3/
    # to the compiler in the system-HepMC3 case, breaking header isolation and
    # potentially causing MARLEY to compile against the wrong headers.
    HEPMC3_INCDIR     := $(INCLUDE_DIR)/builtin
    HEPMC3_CXXFLAGS   := -I$(HEPMC3_INCDIR)
    HEPMC3_LDFLAGS    := -L$(BUILD_LIBDIR) -lMARLEYHepMC3
    HEPMC3_SHARED_LIB := $(BUILD_LIBDIR)/libMARLEYHepMC3.$(SHARED_LIB_SUFFIX)

    # Disable warnings about deprecated declarations in the built-in HepMC3
    # library (triggered via use of sprintf)
    marley_hepmc3.o: HEPMC3_CXXFLAGS += -Wno-deprecated-declarations -Wno-strict-aliasing
  endif

  # The user may force the Makefile to ignore ROOT entirely by defining
  # IGNORE_ROOT="yes" (or any non-empty string) on the command line
  # invocation of make.
  ifndef IGNORE_ROOT
    ROOTCONFIG := $(shell command -v root-config 2> /dev/null)
    # prefer rootcling as the dictionary generator executable name, but use
    # rootcint if you can't find it
    ROOTCLING := $(shell command -v rootcling 2> /dev/null)
    ifndef ROOTCLING
      ROOTCLING := $(shell command -v rootcint 2> /dev/null)
    endif
    ROOT := $(shell command -v root 2> /dev/null)

    ifndef ROOTCONFIG
      $(info WARNING: Could not find a valid ROOT installation.)
      $(info MARLEY will be built without ROOT support.)
      USE_ROOT = no
    else
      ROOT_VERSION := $(shell $(ROOTCONFIG) --version)
      $(info Found ROOT version $(ROOT_VERSION) in $(ROOT))
      $(info MARLEY will be built with ROOT support.)
      override CXXFLAGS += -DUSE_ROOT
      ROOT_CXXFLAGS := $(shell $(ROOTCONFIG) --cflags)

      # The flags MARLEY *publishes* for ROOT, as distinct from the link line
      # MARLEY uses for itself. These are root-config's own output, captured
      # once and emitted verbatim thereafter, so that marley-config and
      # marley.pc describe the ROOT that MARLEY was compiled against rather
      # than whatever root-config says at the moment a consumer runs them.
      #
      # -std= is filtered out because cxx_std is authoritative: it is the
      # maximum of MARLEY's floor, any requested standard, and ROOT's, so
      # ROOT's own -std can only duplicate it or contradict it.
      ROOT_PUBLIC_CFLAGS := $(filter-out -std=%,$(ROOT_CXXFLAGS))
      ROOT_PUBLIC_LIBS := $(shell $(ROOTCONFIG) --libs)

      # The same flags again with ROOT's own -Wl,-rpath removed, for
      # marley.pc. root-config emits an rpath pointing at ROOT's library
      # directory; marley-config may carry it, because a hand-written
      # Makefile has no other way to make its binary runnable. A .pc file
      # must not: a packager's toolchain manages rpath itself, and this one
      # would bake the build machine's ROOT location into a redistributed
      # binary -- working on that machine and nowhere else (spec 6.3).
      ROOT_PC_LIBS := $(filter-out -Wl$(comma)-rpath$(comma)%,$(ROOT_PUBLIC_LIBS))

      ROOT_LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
      ROOT_LIBDIR := $(shell $(ROOTCONFIG) --libdir)
      ROOT_LDFLAGS += -L$(ROOT_LIBDIR) -lCore -lRIO -lHist -lTree
      ifeq ($(UNAME_S),Linux)
        ROOT_LDFLAGS += -rdynamic
      endif
      ROOT_DICT_INCLUDES := -I$(INCLUDE_DIR) -I$(HEPMC3_INCDIR) \
        $(patsubst $(HEPMC3_INCDIR)/%,%,$(wildcard $(HEPMC3_INCDIR)/HepMC3/Data/*.h)) \
        marley/marley_linkdef.hh
      USE_ROOT = yes
      ROOT_OBJ_DICT = marley_root_dict.o

      OBJECTS += marley_root.o OutputFileRoot.o
      OBJECTS += $(ROOT_OBJ_DICT)

$(ROOT_OBJ_DICT):
	$(RM) marley_root_dict*.* libMARLEY_rdict.pcm libMARLEY.rootmap
	# NOTE: "-s <lib>" alone embeds the owning-library name in the .pcm but
	# does not, by itself, guarantee a .rootmap file is written (this varies
	# across ROOT/rootcling versions). "-rml <lib>" and "-rmf <file>" are the
	# same flags ROOT's own ROOT_GENERATE_DICTIONARY() CMake macro uses
	# internally (see cmake/modules/RootMacros.cmake upstream) to force
	# generation of a correctly-named .rootmap, so we mirror that here for
	# robustness across ROOT versions. (Note: this requires a Cling-based
	# "rootcling" from ROOT6+; the legacy ROOT5 "rootcint" fallback used
	# above does not support these flags and is not supported for ROOT
	# dictionary autoloading.)
	$(ROOTCLING) -f marley_root_dict.cc -c \
	  -s $(SHARED_LIB_FILE) -rml $(SHARED_LIB_FILE) -rmf libMARLEY.rootmap \
	  $(ROOT_DICT_INCLUDES)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) \
	  -I$(INCLUDE_DIR) -fPIC \
          -o $(ROOT_OBJ_DICT) -c marley_root_dict.cc
	@mkdir -p $(BUILD_LIBDIR)
	# NOTE: because "-s $(SHARED_LIB_FILE)" is passed above, rootcling names
	# the generated pcm/rootmap files after that library (libMARLEY), not
	# after the "marley_root_dict" dictionary source file name.
	mv $(BUILD_DIR)/libMARLEY_rdict.pcm $(BUILD_LIBDIR)/ 2>/dev/null || true
	mv $(BUILD_DIR)/libMARLEY.rootmap $(BUILD_LIBDIR)/ 2>/dev/null || true
	$(RM) marley_root_dict.cc

    endif
  else
    $(info Ignoring any ROOT installations that may be present.)
    $(info MARLEY will be built without ROOT support.)
    USE_ROOT = no
  endif
endif
endif

# ---------------------------------------------------------------------------
# The effective C++ standard (design spec 6.2)
#
# One value: the maximum of MARLEY's floor of 17, whatever the user asked
# for, and whatever ROOT was built with. It is used for compilation and is
# what marley-config --cxx-std and marley.pc report, so those cannot
# disagree with the binary.
#
# This replaces four findstring blocks that mis-parsed root-config's output
# in two separate ways. They could DOWNGRADE -- "make CXX_STD=c++23" against
# a C++20 ROOT fired the c++20 block and overrode the request back down --
# and they were blind to the gnu++ spellings, since "-std=gnu++20" does not
# contain the substring "c++20" and root-config emits gnu++ on at least one
# platform branch.
# ---------------------------------------------------------------------------

# gnu++NN and c++NN alike reduce to NN, and the provisional spellings map to
# their final numbers. The gnu++ substitution runs first, so that "gnu++2a"
# reaches the c++2a rule at all.
marley_std_num = $(patsubst c++%,%,$(subst c++2c,c++26,$(subst c++2b,c++23,\
  $(subst c++2a,c++20,$(subst gnu++,c++,$(1))))))

# Sorted ascending, so the largest surviving entry is the last word.
STD_ORDER := 17 20 23 26

ROOT_STD_RAW := $(patsubst -std=%,%,\
  $(lastword $(filter -std=%,$(ROOT_CXXFLAGS))))

CXX_STD_REQ_NUM := $(call marley_std_num,$(CXX_STD))
ROOT_STD_NUM := $(call marley_std_num,$(ROOT_STD_RAW))

# filter keeps the words of STD_ORDER that appear among the candidates, so an
# unrecognised standard drops out rather than winning; 17 is always a
# candidate, so the result can never fall below the floor.
CXX_STD_CANDIDATES := 17 $(CXX_STD_REQ_NUM) $(ROOT_STD_NUM)
CXX_STD_NUM := $(lastword $(filter $(CXX_STD_CANDIDATES),$(STD_ORDER)))

CXX_STD_BEFORE := $(CXX_STD)
override CXX_STD := c++$(CXX_STD_NUM)

ifneq ($(CXX_STD),$(CXX_STD_BEFORE))
  ifdef CXX_STD_EXPLICIT
    # The user named a standard and is not getting it. That is worth a
    # warning rather than a status line.
    $(warning Building with -std=$(CXX_STD), not the requested $(CXX_STD_BEFORE)$(if $(filter $(CXX_STD_NUM),$(ROOT_STD_NUM)), -- ROOT was built with it, -- it is MARLEY's minimum))
  else
    $(info Using -std=$(CXX_STD) to match ROOT.)
  endif
  # Appended rather than edited in: CXXFLAGS was assembled before ROOT was
  # found, and a later -std takes precedence over an earlier one.
  override CXXFLAGS += -std=$(CXX_STD)
endif

# (The install directory variables are defined near the top of this file,
# because the build tree mirrors the install layout and SHARED_LIB needs
# libdir_rel.)

MARLEY_LIBS := $(SHARED_LIB)
ifeq ($(USE_ROOT),yes)
  # If we're building the tests, then link to some extra ROOT libraries
  ifneq (,$(findstring test,$(MAKECMDGOALS)))
    ROOT_LDFLAGS += -lGraf -lGpad
  endif
endif

# Causes GNU make to auto-delete the object files when the build is complete
.INTERMEDIATE: $(OBJECTS) $(TEST_OBJECTS) marley_hepmc3.o marley.o

# Define a variable to store a command used repeatedly in target definitions
COMPILE_CXX = $(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) \
              $(HEPMC3_CXXFLAGS) -I$(INCLUDE_DIR) -fPIC

%.o: $(SRC_DIR)/%.cc
	$(COMPILE_CXX) -o $@ -c $<

%.o: $(SRC_DIR)/app/%.cc
	$(COMPILE_CXX) -o $@ -c $<

%.o: $(SRC_DIR)/tests/%.cc
	$(COMPILE_CXX) -o $@ -c $<

# Special rule used as needed to compile the built-in GSL code without
# hidden symbols so we can use them in the test executable
%.test.o: $(SRC_DIR)/%.cc
	$(COMPILE_CXX) -o $@ -c $<


ifndef FOUND_HEPMC3

marley_hepmc3.o: $(SRC_DIR)/marley_hepmc3.cc
	$(COMPILE_CXX) -o $@ -c $<

$(HEPMC3_SHARED_LIB): marley_hepmc3.o
	@mkdir -p $(BUILD_LIBDIR)
	$(CXX) $(CXXFLAGS) -fPIC -shared -o $@ $< \
	  $(HEPMC3_INSTALL_NAME_FLAG) $(NEW_DTAGS) -lm

endif

# Make the build tree mirror the install layout, so that the single hop taken
# by self-location (see the design spec, 2.2) works before installation and
# the whole "works in the build tree, breaks when installed" class of defects
# cannot arise. Absolute link targets; -n so that re-running does not create a
# link *inside* an existing one. An order-only prerequisite rather than a
# one-time setup step, because `clean` is `$(RM) -rf $(BUILD_DIR)`.
BUILD_TREE_LINKS := $(BUILD_DIR)/include $(BUILD_DIR)/share/marley/data

$(BUILD_DIR)/include:
	@mkdir -p $(BUILD_DIR)
	@ln -sfn $(abspath $(INCLUDE_DIR)) $@

$(BUILD_DIR)/share/marley/data:
	@mkdir -p $(BUILD_DIR)/share/marley
	@ln -sfn $(abspath $(TOP_DIR)/data) $@

$(SHARED_LIB): $(HEPMC3_SHARED_LIB) $(OBJECTS) | $(BUILD_TREE_LINKS)
	@mkdir -p $(BUILD_LIBDIR)
	$(CXX) $(CXXFLAGS) $(ROOT_CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) \
	  -fPIC -shared -o $@ $(OBJECTS) \
	  $(GSL_LDFLAGS) $(ROOT_LDFLAGS) $(HEPMC3_LDFLAGS) $(DL_LIBS) \
	  $(INSTALL_NAME_FLAG) $(NEW_DTAGS) $(RPATH_LIB) $(EXT_RPATH)

marley: $(BUILD_DIR)/bin/marley

mroot: $(BUILD_DIR)/bin/mroot

setup-script: $(BUILD_DIR)/bin/setup_marley.sh

marley-config: $(BUILD_DIR)/bin/marley-config

# ---------------------------------------------------------------------------
# Generated shell scripts (design spec 6.4)
#
# mroot and setup_marley.sh are each generated TWICE from one template: once
# for the build tree and once for the install tree. The two trees differ by
# more than a prefix -- a vendored HepMC3 sits at include/builtin before
# installation and at include/marley/builtin after it -- and a single
# self-locating script applying one set of relative hops cannot express that.
#
# The two copies share a basename, so the install-tree copies are staged in
# their own directory and installed from there.
# ---------------------------------------------------------------------------
INSTALL_SCRIPT_DIR := $(BUILD_DIR)/install-scripts

# The build tree's include symlink is named "include" unconditionally, while
# the install tree honours includedir. Likewise bin.
BUILD_BIN_REL     := bin
BUILD_INCLUDE_REL := include
BUILD_BUILTIN_REL  = $(if $(FOUND_HEPMC3),,include/builtin)
INSTALL_BUILTIN_REL = $(if $(FOUND_HEPMC3),,$(includedir_rel)/marley/builtin)

# $(1) bin_rel  $(2) include_rel  $(3) builtin_rel
define gen_script
	@mkdir -p $(dir $@)
	sed -e 's|@bin_rel@|$(1)|g' \
	  -e 's|@lib_rel@|$(libdir_rel)|g' \
	  -e 's|@include_rel@|$(2)|g' \
	  -e 's|@builtin_rel@|$(3)|g' \
	  -e 's|@hepmc3_include@|$(if $(FOUND_HEPMC3),$(HEPMC3_INCDIR))|g' \
	  $< > $@
	chmod +x $@
endef

$(BUILD_DIR)/bin/mroot: $(SRC_DIR)/scripts/mroot.in
	$(call gen_script,$(BUILD_BIN_REL),$(BUILD_INCLUDE_REL),$(BUILD_BUILTIN_REL))

$(BUILD_DIR)/bin/setup_marley.sh: $(SRC_DIR)/scripts/setup_marley.sh.in
	$(call gen_script,$(BUILD_BIN_REL),$(BUILD_INCLUDE_REL),$(BUILD_BUILTIN_REL))

$(INSTALL_SCRIPT_DIR)/mroot: $(SRC_DIR)/scripts/mroot.in
	$(call gen_script,$(bindir_rel),$(includedir_rel),$(INSTALL_BUILTIN_REL))

$(INSTALL_SCRIPT_DIR)/setup_marley.sh: $(SRC_DIR)/scripts/setup_marley.sh.in
	$(call gen_script,$(bindir_rel),$(includedir_rel),$(INSTALL_BUILTIN_REL))

# marley-config is concatenated, not substituted (design spec 6.2): a
# generated header of KEY='value' assignments, followed by the static body
# shipped in src/scripts/. The old mechanism transcribed an eight-token
# substitution list into two languages and failed silently when a token was
# missed, because an unsubstituted @@ROOT_LIBS@@ is not a shell syntax
# error. There is now no mapping to get wrong, and the body's :? contract
# catches an incomplete header in one place.
#
# Note what is NOT here: no MARLEY-owned path. Those are computed by the
# body from the script's own location, so this rule is structurally unable
# to bake in a build-tree -L or a source-tree -I.
$(BUILD_DIR)/bin/marley-config: $(SRC_DIR)/scripts/marley-config-body.sh $(MARLEY_LIBS)
	@mkdir -p $(BUILD_DIR)/bin
	printf '%s\n' \
	  '#!/bin/bash' \
	  '# Generated at build time by make/build.mk -- do not edit.' \
	  '# The static half of this script is src/scripts/marley-config-body.sh.' \
	  "version='$(MARLEY_VERSION)'" \
	  "git_revision='$(GIT_REVISION)'" \
	  "cxx_std='$(CXX_STD)'" \
	  "compiled_log_level='$(MARLEY_LOG_LEVEL)'" \
	  "use_root='$(USE_ROOT)'" \
	  "hepmc3_mode='$(if $(FOUND_HEPMC3),external,builtin)'" \
	  "gsl_mode='$(if $(FOUND_GSL),external,builtin)'" \
	  "root_cflags='$(ROOT_PUBLIC_CFLAGS)'" \
	  "root_libs='$(ROOT_PUBLIC_LIBS)'" \
	  "ext_hepmc3_cflags='$(HEPMC3_PUBLIC_CFLAGS)'" \
	  "ext_hepmc3_libs='$(HEPMC3_PUBLIC_LIBS)'" \
	  "ext_gsl_cflags='$(if $(FOUND_GSL),$(GSL_CXXFLAGS))'" \
	  "ext_gsl_libs='$(if $(FOUND_GSL),$(GSL_LDFLAGS))'" \
	  "new_dtags='$(NEW_DTAGS)'" \
	  "libdir_rel='$(libdir_rel)'" \
	  "bindir_rel='$(bindir_rel)'" \
	  "includedir_rel='$(includedir_rel)'" \
	  > $@
	cat $(SRC_DIR)/scripts/marley-config-body.sh >> $@
	chmod +x $@

# marley.pc self-locates from ${pcfiledir}, so like marley-config it needs
# only the relative directory names. The conditional fields are whole-string
# substitution keys rather than any kind of bracket notation, because
# configure_file (the CMake side) cannot express a conditional and both
# sides must produce the same file.
#
# GSL appears only as Requires.private: it is used by src/marley_gsl.cc and
# named in no installed header, so a consumer neither needs nor sees it.
# HepMC3 is public -- MARLEY's event type is an HepMC3 type.
# marley.pc carries the same flags marley-config does, minus the rpath. It
# has no Requires: line, because neither HepMC3 nor ROOT ships a pkg-config
# file and GSL is private -- see the template for the full reasoning.
#
# -DUSE_ROOT is not separable from ROOT's include path: EventFileReader.hh
# includes TFile.h and TTree.h under that macro, so publishing the macro
# without the -I hands the consumer a public header that cannot be parsed.
PC_CFLAGS_EXTRA := $(strip \
  $(if $(filter yes,$(USE_ROOT)),-DUSE_ROOT $(ROOT_PUBLIC_CFLAGS)) \
  $(if $(FOUND_HEPMC3),$(HEPMC3_PUBLIC_CFLAGS),-I$${includedir}/marley/builtin))
PC_LIBS_EXTRA := $(strip \
  $(if $(FOUND_HEPMC3),$(HEPMC3_PUBLIC_LIBS),-lMARLEYHepMC3) \
  $(if $(filter yes,$(USE_ROOT)),$(ROOT_PC_LIBS)))

marley-pc: $(BUILD_LIBDIR)/pkgconfig/marley.pc

$(BUILD_LIBDIR)/pkgconfig/marley.pc: $(SRC_DIR)/scripts/marley.pc.in $(MARLEY_LIBS)
	@mkdir -p $(dir $@)
	sed -e 's|@libdir_rel@|$(libdir_rel)|g' \
	  -e 's|@includedir_rel@|$(includedir_rel)|g' \
	  -e 's|@version@|$(MARLEY_VERSION)|g' \
	  -e 's|@cxx_std@|$(CXX_STD)|g' \
	  -e 's|@pc_cflags_extra@|$(PC_CFLAGS_EXTRA)|g' \
	  -e 's|@pc_libs_extra@|$(PC_LIBS_EXTRA)|g' \
	  $< > $@

$(BUILD_DIR)/bin/marley: $(MARLEY_LIBS) marley.o $(BUILD_DIR)/bin/marley-config $(if $(filter yes,$(USE_ROOT)),$(BUILD_DIR)/bin/mroot)
	@mkdir -p $(BUILD_DIR)/bin
	$(CXX) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) -o $@ marley.o \
	  -L$(BUILD_LIBDIR) -l$(SHARED_LIB_NAME) \
	  $(ROOT_LDFLAGS) $(if $(filter yes,$(USE_ROOT)),-lGraf -lGpad) \
	  $(GSL_LDFLAGS) $(HEPMC3_LDFLAGS) $(NEW_DTAGS) $(RPATH_EXE) $(EXT_RPATH)

$(TEST_EXECUTABLE): $(TEST_OBJECTS) $(MARLEY_LIBS)
	@mkdir -p $(BUILD_DIR)/bin
	$(CXX) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) -o $@ $(TEST_OBJECTS) \
	  -L$(BUILD_LIBDIR) -l$(SHARED_LIB_NAME) $(ROOT_LDFLAGS) \
	  $(GSL_LDFLAGS) $(HEPMC3_LDFLAGS) $(NEW_DTAGS) $(RPATH_EXE) $(EXT_RPATH)

# Build the marg4 Geant4 example (requires a separate Geant4 installation).
# Delegates to the example's own Makefile with PATH set so that it finds this
# build's marley-config, which then reports this build's directories by
# locating itself. No MARLEY variable is exported: the example no longer has
# one to read, and that is the point of it being an example.
#
# BINDIR/OBJDIR are passed so that an in-tree build still lands in
# build/bin alongside marley; an installed copy of the example defaults to
# building beside itself instead.
.PHONY: marg4
marg4: $(SHARED_LIB) $(BUILD_DIR)/bin/marley-config
	@PATH=$(BUILD_DIR)/bin:$$PATH \
	  $(MAKE) -C $(TOP_DIR)/examples/marg4 \
	    BINDIR=$(BUILD_DIR)/bin OBJDIR=$(BUILD_DIR)/marg4_obj

DOCS_BUILD_DIR = $(BUILD_DIR)/docs

.PHONY: marley mroot setup-script marley-config marley-pc docs clean install uninstall

doxygen:
	export MARLEY_VERSION=$(VERSION_PREFIX)$(MARLEY_VERSION) \
	export DOXYGEN_OUTPUT_DIR=$(DOCS_BUILD_DIR)/doxygen \
	export BUILDDIR=$(DOCS_BUILD_DIR) \
	&& cd $(TOP_DIR)/docs && $(MAKE) doxygen

docs:
	export MARLEY_VERSION=$(VERSION_PREFIX)$(MARLEY_VERSION) \
	export DOXYGEN_OUTPUT_DIR=$(DOCS_BUILD_DIR)/doxygen \
	export BUILDDIR=$(DOCS_BUILD_DIR) \
	&& cd $(TOP_DIR)/docs && $(MAKE) html

clean:
	$(RM) -rf $(BUILD_DIR)

install: marley marley-config marley-pc \
  $(INSTALL_SCRIPT_DIR)/setup_marley.sh \
  $(if $(filter yes,$(USE_ROOT)),$(INSTALL_SCRIPT_DIR)/mroot)
	mkdir -p $(DESTDIR)$(bindir)
	mkdir -p $(DESTDIR)$(libdir)
	mkdir -p $(DESTDIR)$(includedir)/marley
	mkdir -p $(DESTDIR)$(datadir)/marley
	cp $(BUILD_DIR)/bin/marley $(DESTDIR)$(bindir)
	cp $(BUILD_DIR)/bin/marley-config $(DESTDIR)$(bindir)
	# From the staging directory, not build/bin: the two copies share a
	# basename and differ in the hops they take (design spec 6.4).
	if [ "$(USE_ROOT)" = "yes" ]; then \
	  cp $(INSTALL_SCRIPT_DIR)/mroot $(DESTDIR)$(bindir); fi
	cp $(INSTALL_SCRIPT_DIR)/setup_marley.sh $(DESTDIR)$(bindir)
	cp $(SHARED_LIB) $(DESTDIR)$(libdir)
	cp $(BUILD_LIBDIR)/libMARLEY_rdict.pcm $(DESTDIR)$(libdir) 2> /dev/null || true
	cp $(BUILD_LIBDIR)/libMARLEY.rootmap $(DESTDIR)$(libdir) 2> /dev/null || true
	# pkgconfig follows libdir: on a lib64 build this is lib64/pkgconfig, and
	# the file's own ${pcfiledir}/../.. still lands on the prefix because
	# libdir_rel is one segment by construction.
	mkdir -p $(DESTDIR)$(libdir)/pkgconfig
	cp $(BUILD_LIBDIR)/pkgconfig/marley.pc $(DESTDIR)$(libdir)/pkgconfig
	cp -r $(TOP_DIR)/data $(DESTDIR)$(datadir)/marley
	# data/tests/ is used only by martest, which runs from the build tree
	rm -rf $(DESTDIR)$(datadir)/marley/data/tests
	cp -r $(TOP_DIR)/include/marley $(DESTDIR)$(includedir)
	# examples/ ships: the documentation references it 15 times across 5
	# files, including commands users are told to type verbatim, and
	# examples/config/annotated.js is the primary job-configuration
	# reference. None of it worked on an installed MARLEY.
	cp -r $(TOP_DIR)/examples $(DESTDIR)$(datadir)/marley
	# COPYING and GUIDELINES ship because the licence header of *every
	# installed header* points at them, and MARLEY is GPLv3. CITATION.bib
	# joins them because the documentation asks users to cite.
	cp $(TOP_DIR)/COPYING $(TOP_DIR)/GUIDELINES $(TOP_DIR)/CITATION.bib \
	  $(DESTDIR)$(datadir)/marley
ifndef FOUND_HEPMC3
	# The vendored HepMC3 is installed under MARLEY's own names: a
	# libHepMC3 of ours in a shared prefix collides with a real one, and
	# this amalgamated copy does not implement all of HepMC3. Nesting the
	# headers changes no #include directive -- the -I root travels with
	# them (design spec 4.1).
	cp $(HEPMC3_SHARED_LIB) $(DESTDIR)$(libdir)
	mkdir -p $(DESTDIR)$(includedir)/marley/builtin/HepMC3
	cp -r $(INCLUDE_DIR)/builtin/HepMC3/. \
	  $(DESTDIR)$(includedir)/marley/builtin/HepMC3
endif
	# The install receipt (design spec 7.1): a coarse manifest recorded in
	# the thing it describes, so it cannot go stale and a CMake install is
	# covered by the same reader. Paths are prefix-relative, so a moved
	# installation still uninstalls. It replaces a 60-line manifest writer
	# and the .marley-installed-builtin-hepmc3 sentinel, which recorded the
	# same fact somewhere else.
	#
	# The glob is libMARLEY*, not libMARLEY.* -- the latter does not match
	# libMARLEY_rdict.pcm, which has an underscore, so a ROOT-enabled
	# install would quietly leave the dictionary behind.
	@mkdir -p $(DESTDIR)$(datadir)/marley
	{ \
	  printf '%s\n' "marley-install-receipt: 1"; \
	  printf '%s\n' "version: $(MARLEY_VERSION)"; \
	  printf '%s\n' "git-revision: $(GIT_REVISION)"; \
	  printf '%s\n' "prefix: $(prefix)"; \
	  printf '%s\n' "build-system: make"; \
	  printf '%s\n' "root: $(USE_ROOT)"; \
	  printf '%s\n' "hepmc3: $(if $(FOUND_HEPMC3),external,builtin)"; \
	  printf '%s\n' "gsl: $(if $(FOUND_GSL),external,builtin)"; \
	  printf '%s\n' "trees: $(includedir_rel)/marley share/marley"; \
	  printf '%s\n' "files: $(bindir_rel)/marley $(bindir_rel)/marley-config $(bindir_rel)/setup_marley.sh"; \
	  if [ "$(USE_ROOT)" = "yes" ]; then \
	    printf '%s\n' "       $(bindir_rel)/mroot"; \
	  fi; \
	  printf '%s\n' "       $(libdir_rel)/libMARLEY*"; \
	  printf '%s\n' "       $(libdir_rel)/pkgconfig/marley.pc"; \
	} > $(DESTDIR)$(datadir)/marley/install-receipt
	$(LDCONFIG_CMD)

# Delegates to the same script the CMake path and a direct invocation use.
# One implementation, three entry points (design spec 7.2).
uninstall:
	$(SRC_DIR)/scripts/uninstall-marley --prefix "$(DESTDIR)$(prefix)"

