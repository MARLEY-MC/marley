# Names of various MARLEY directories
TOP_DIR ?= $(CURDIR)
BUILD_DIR = $(TOP_DIR)/build
DATA_DIR = $(TOP_DIR)/data
INCLUDE_DIR = $(TOP_DIR)/include
SRC_DIR = $(TOP_DIR)/src

# Default version of the C++ standard to use for compilation
CXX_STD=c++17

UNKNOWN_REV=unknown version

# If the target is "debug", then include debugging information and build the
# marley executable with optimization turned off. Prepend the default -std
# option so that it will be overridden if the user has manually specified one
# in CXXFLAGS.
ifeq ($(MAKECMDGOALS),debug)
  override CXXFLAGS := -std=$(CXX_STD) $(CXXFLAGS) -O0 -g \
    -DMARLEY_COMPILED_LOG_LEVEL=TRACE
else
  # Otherwise, use full optimization and do not include debugging info. Prepend
  # the default -O3 and -std options here in case the user wants to manually
  # override them in CXXFLAGS.
  override CXXFLAGS := -O3 -std=$(CXX_STD) $(CXXFLAGS) \
    -DMARLEY_COMPILED_LOG_LEVEL=INFO
endif

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SHARED_LIB_SUFFIX=dylib
else ifeq ($(UNAME_S),Linux)
  SHARED_LIB_SUFFIX=so
else
  $(warning Unrecognized operating system encountered.)
  SHARED_LIB_SUFFIX=so
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

  # Define a link to a tarball on GitHub for the current tagged release
  TARBALL_LINK = "<a href=\"https://github.com/MARLEY-MC/marley/$\
    archive/v$(MARLEY_VERSION).tar.gz\">here</a> or"

  VERSION_PREFIX="v"

else

  # Just link to the repository in the doxygen documentation if we're not
  # working with a tagged release
  TARBALL_LINK = ""

  # Also use the git revision as the version number in this case
  MARLEY_VERSION := $(GIT_REVISION)

endif

# Define the MARLEY_VERSION preprocessor macro
override CXXFLAGS += -DMARLEY_VERSION="\"$(MARLEY_VERSION)\""

# Define the MARLEY_GIT_REVISION preprocessor macro
override CXXFLAGS += -DMARLEY_GIT_REVISION="\"$(GIT_REVISION)\""

SHARED_LIB_NAME := MARLEY
SHARED_LIB_FILE := lib$(SHARED_LIB_NAME).$(SHARED_LIB_SUFFIX)
SHARED_LIB := $(BUILD_DIR)/lib/$(SHARED_LIB_FILE)

TEST_EXECUTABLE = $(BUILD_DIR)/bin/martest
TEST_OBJECTS = $(notdir $(patsubst %.cc,%.o,$(wildcard $(SRC_DIR)/tests/*.cc)))

all: marley
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
      ifneq (,$(findstring $(COMPILER_VERSION), 4.9.))
        # g++ 4.9 gives many false positives for -Wshadow, so disable it
        # for now.
        override CXXFLAGS += -Wno-shadow
      endif
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
  OBJECTS := $(filter-out OutputFilePlainRoot.o, $(OBJECTS))
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
      FOUND_HEPMC3    := TRUE
    endif
  else
    $(info Ignoring any HepMC3 installations that may be present.)
  endif

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
    HEPMC3_LDFLAGS    := -L$(BUILD_DIR)/lib -lHepMC3
    HEPMC3_SHARED_LIB := $(BUILD_DIR)/lib/libHepMC3.$(SHARED_LIB_SUFFIX)

    # Disable warnings about deprecated declarations in the built-in HepMC3
    # library (triggered via use of sprintf)
    marley_hepmc3.o: HEPMC3_CXXFLAGS += -Wno-deprecated-declarations
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
      CXXFLAGS += -std=$(CXX_STD)
    else
      ROOT_VERSION := $(shell $(ROOTCONFIG) --version)
      $(info Found ROOT version $(ROOT_VERSION) in $(ROOT))
      $(info MARLEY will be built with ROOT support.)
      override CXXFLAGS += -DUSE_ROOT
      ROOT_CXXFLAGS := $(shell $(ROOTCONFIG) --cflags)

      # If ROOT was built with a later C++ standard, switch to building MARLEY
      # with it as well, just in case. Later -std options take precedence, so
      # just tack on a new one to keep things simple.
      ifneq (, $(findstring c++2a, $(ROOT_CXXFLAGS)))
        override CXX_STD = c++2a
        override CXXFLAGS += -std=$(CXX_STD)
      endif
      ifneq (, $(findstring c++20, $(ROOT_CXXFLAGS)))
        override CXX_STD = c++20
        override CXXFLAGS += -std=$(CXX_STD)
      endif
      ifneq (, $(findstring c++2b, $(ROOT_CXXFLAGS)))
        override CXX_STD = c++2b
        override CXXFLAGS += -std=$(CXX_STD)
      endif
      ifneq (, $(findstring c++23, $(ROOT_CXXFLAGS)))
        override CXX_STD = c++23
        override CXXFLAGS += -std=$(CXX_STD)
      endif

      ROOT_LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
      ROOT_LIBDIR := $(shell $(ROOTCONFIG) --libdir)
      ROOT_LDFLAGS += -L$(ROOT_LIBDIR) -lCore -lRIO -lHist -lTree
      ifeq ($(UNAME_S),Linux)
        ROOT_LDFLAGS += -rdynamic
      endif
      ROOT_DICT_INCLUDES := -I$(INCLUDE_DIR) -I$(HEPMC3_INCDIR) \
        -I$(HEPMC3_INCDIR)/HepMC3/Data \
        $(notdir $(wildcard $(HEPMC3_INCDIR)/HepMC3/Data/*.h)) \
        marley/marley_linkdef.hh
      USE_ROOT = yes
      ROOT_OBJ_DICT = marley_root_dict.o

      OBJECTS += marley_root.o OutputFileRoot.o
      OBJECTS += OutputFilePlainRoot.o $(ROOT_OBJ_DICT)

$(ROOT_OBJ_DICT):
	$(RM) marley_root_dict*.*
	$(ROOTCLING) -f marley_root_dict.cc -c $(ROOT_DICT_INCLUDES)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) \
	  -I$(INCLUDE_DIR) -I$(HEPMC3_INCDIR)/HepMC3/Data -fPIC \
          -o $(ROOT_OBJ_DICT) -c marley_root_dict.cc
	@mkdir -p $(BUILD_DIR)/lib
	mv $(BUILD_DIR)/marley_root_dict_rdict.pcm $(BUILD_DIR)/lib/ 2>/dev/null || true
	mv $(BUILD_DIR)/marley_root_dict.rootmap $(BUILD_DIR)/lib/ 2>/dev/null || true
	$(RM) marley_root_dict.cc

    endif
  else
    $(info Ignoring any ROOT installations that may be present.)
    $(info MARLEY will be built without ROOT support.)
    USE_ROOT = no
  endif
endif
endif

# Set up default variables for the install/uninstall targets and for setting
# the executable rpath
prefix = /usr
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
datadir = $(prefix)/share
libdir = $(exec_prefix)/lib
incdir = $(prefix)/include

# Path (without DESTDIR) to the install manifest that records which files were
# placed on disk by "make install".  This is the single source of truth used by
# "make uninstall" to remove exactly the installed files.
MARLEY_MANIFEST_RELPATH = $(libdir)/marley-install-manifest.txt
MARLEY_MANIFEST = $(DESTDIR)$(MARLEY_MANIFEST_RELPATH)

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
	@mkdir -p $(BUILD_DIR)/lib
	$(CXX) $(CXXFLAGS) -fPIC -shared -o $@ $<

endif

$(SHARED_LIB): $(HEPMC3_SHARED_LIB) $(OBJECTS)
	@mkdir -p $(BUILD_DIR)/lib
	$(CXX) $(CXXFLAGS) $(ROOT_CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) $(GSL_LDFLAGS) $(ROOT_LDFLAGS) $(HEPMC3_LDFLAGS) \
	-fPIC -shared -o $@ $(OBJECTS)

marley: $(BUILD_DIR)/bin/marley

mroot: $(BUILD_DIR)/bin/mroot

marley-config: $(BUILD_DIR)/bin/marley-config

$(BUILD_DIR)/bin/mroot: $(MARLEY_LIBS)
	@mkdir -p $(BUILD_DIR)/bin
	cp $(SRC_DIR)/scripts/mroot $@

# We use a temporary backup file here so that the invocation of sed is
# compatible with both the GNU/Linux and BSD/macOS versions.
# See https://stackoverflow.com/a/22084103/4081973 for details.
$(BUILD_DIR)/bin/marley-config: $(MARLEY_LIBS)
	@mkdir -p $(BUILD_DIR)/bin
	$(RM) $(BUILD_DIR)/bin/marley-config
	cp $(SRC_DIR)/scripts/marley-config.in $(BUILD_DIR)/bin/marley-config
	sed -i.bak -e '/^##/d' -e "s|@@VERSION@@|\"$(MARLEY_VERSION)\"|g" \
	  -e "s|@@GIT_REVISION@@|\"$(GIT_REVISION)\"|g" \
	  -e "s|@@CXX_STD@@|\"$(CXX_STD)\"|g" \
	  -e "s|@@HEPMC3_CFLAGS@@|\"$(HEPMC3_CXXFLAGS)\"|g" \
	  -e "s|@@HEPMC3_LIBS@@|\"$(HEPMC3_LDFLAGS)\"|g" \
	  -e "s|@@GSL_CFLAGS@@|\"$(GSL_CXXFLAGS)\"|g" \
	  -e "s|@@GSL_LIBS@@|\"$(GSL_LDFLAGS)\"|g" \
	  -e "s|@@USE_ROOT@@|\"$(USE_ROOT)\"|g" $(BUILD_DIR)/bin/marley-config
	$(RM) $(BUILD_DIR)/bin/marley-config.bak

$(BUILD_DIR)/bin/marley: $(MARLEY_LIBS) marley.o $(BUILD_DIR)/bin/marley-config $(if $(filter yes,$(USE_ROOT)),$(BUILD_DIR)/bin/mroot)
	@mkdir -p $(BUILD_DIR)/bin
	$(CXX) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) -o $@ -L$(BUILD_DIR)/lib \
	  -l$(SHARED_LIB_NAME) $(ROOT_LDFLAGS) $(if $(filter yes,$(USE_ROOT)),-lGraf -lGpad) \
	  $(GSL_LDFLAGS) $(HEPMC3_LDFLAGS) -Wl,-rpath -Wl,$(libdir):$(BUILD_DIR)/lib marley.o

$(TEST_EXECUTABLE): $(TEST_OBJECTS) $(MARLEY_LIBS)
	@mkdir -p $(BUILD_DIR)/bin
	$(CXX) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) -o $@ -L$(BUILD_DIR)/lib \
	  -l$(SHARED_LIB_NAME) $(ROOT_LDFLAGS) \
	  $(GSL_LDFLAGS) $(HEPMC3_LDFLAGS) $(TEST_OBJECTS)

.PHONY: marley mroot docs clean install uninstall

doxygen:
	export MARLEY_VERSION=$(VERSION_PREFIX)$(MARLEY_VERSION) \
	export TARBALL_LINK=$(TARBALL_LINK) \
	&& cd $(TOP_DIR)/docs && $(MAKE) doxygen

docs:
	export MARLEY_VERSION=$(VERSION_PREFIX)$(MARLEY_VERSION) \
	export TARBALL_LINK=$(TARBALL_LINK) \
	&& cd $(TOP_DIR)/docs && $(MAKE) html

clean:
	$(RM) -rf $(BUILD_DIR)

install: marley marley-config $(if $(filter yes,$(USE_ROOT)),mroot)
	mkdir -p $(DESTDIR)$(bindir)
	mkdir -p $(DESTDIR)$(libdir)
	mkdir -p $(DESTDIR)$(incdir)/marley
	mkdir -p $(DESTDIR)$(datadir)/marley
	cp $(BUILD_DIR)/bin/marley $(DESTDIR)$(bindir)
	cp $(BUILD_DIR)/bin/marley-config $(DESTDIR)$(bindir)
	if [ "$(USE_ROOT)" = "yes" ]; then cp $(BUILD_DIR)/bin/mroot $(DESTDIR)$(bindir); fi
	cp $(SHARED_LIB) $(DESTDIR)$(libdir)
	cp $(BUILD_DIR)/lib/marley_root_dict_rdict.pcm $(DESTDIR)$(libdir) 2> /dev/null || true
	cp $(BUILD_DIR)/lib/marley_root_dict.rootmap $(DESTDIR)$(libdir) 2> /dev/null || true
	cp -r $(TOP_DIR)/data $(DESTDIR)$(datadir)/marley
	cp -r $(TOP_DIR)/include/marley $(DESTDIR)$(incdir)
ifndef FOUND_HEPMC3
	cp $(HEPMC3_SHARED_LIB) $(DESTDIR)$(libdir)
	mkdir -p $(DESTDIR)$(incdir)/HepMC3
	cp -r $(INCLUDE_DIR)/builtin/HepMC3/. $(DESTDIR)$(incdir)/HepMC3
	touch $(DESTDIR)$(libdir)/.marley-installed-builtin-hepmc3
endif
	# Write the install manifest: one absolute path per line (no DESTDIR prefix),
	# matching the convention used by CMake's install_manifest.txt.  The manifest
	# records every file actually placed on disk so that "make uninstall" can
	# remove them precisely without touching unrelated files or directories.
	{ \
	  printf '%s\n' "$(bindir)/marley"; \
	  printf '%s\n' "$(bindir)/marley-config"; \
	  if [ "$(USE_ROOT)" = "yes" ]; then printf '%s\n' "$(bindir)/mroot"; fi; \
	  printf '%s\n' "$(libdir)/$(SHARED_LIB_FILE)"; \
	  if [ -f "$(DESTDIR)$(libdir)/marley_root_dict_rdict.pcm" ]; then \
	    printf '%s\n' "$(libdir)/marley_root_dict_rdict.pcm"; \
	  fi; \
	  if [ -f "$(DESTDIR)$(libdir)/marley_root_dict.rootmap" ]; then \
	    printf '%s\n' "$(libdir)/marley_root_dict.rootmap"; \
	  fi; \
	  find "$(DESTDIR)$(datadir)/marley" -type f \
	    | sed 's|^$(DESTDIR)||'; \
	  find "$(DESTDIR)$(incdir)/marley" -type f \
	    | sed 's|^$(DESTDIR)||'; \
	  if [ -f "$(DESTDIR)$(libdir)/.marley-installed-builtin-hepmc3" ]; then \
	    printf '%s\n' "$(libdir)/libHepMC3.$(SHARED_LIB_SUFFIX)"; \
	    find "$(DESTDIR)$(incdir)/HepMC3" -type f \
	      | sed 's|^$(DESTDIR)||'; \
	    printf '%s\n' "$(libdir)/.marley-installed-builtin-hepmc3"; \
	  fi; \
	} > $(MARLEY_MANIFEST)
	ldconfig

uninstall:
	@if [ ! -f "$(MARLEY_MANIFEST)" ]; then \
	  echo "ERROR: Install manifest not found: $(MARLEY_MANIFEST)"; \
	  echo "       Have you run 'make install' from this build? Aborting."; \
	  exit 1; \
	fi
	# Remove every file listed in the manifest.
	while IFS= read -r f; do \
	  if [ -f "$(DESTDIR)$$f" ] || [ -L "$(DESTDIR)$$f" ]; then \
	    echo "Removing: $(DESTDIR)$$f"; \
	    rm -f "$(DESTDIR)$$f"; \
	  fi; \
	done < $(MARLEY_MANIFEST)
	# Remove the manifest itself.
	rm -f $(MARLEY_MANIFEST)
	# Prune directories that may now be empty.  Use plain rmdir (POSIX) rather
	# than the GNU-specific --ignore-fail-on-non-empty flag so this works on
	# non-Linux systems (e.g. macOS/BSD).  The '|| true' suppresses the exit
	# code that rmdir emits when a directory is not empty.
	rmdir -p "$(DESTDIR)$(datadir)/marley" 2>/dev/null || true
	rmdir "$(DESTDIR)$(incdir)/marley" 2>/dev/null || true
	rmdir "$(DESTDIR)$(incdir)/HepMC3" 2>/dev/null || true
	ldconfig
