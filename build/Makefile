INCLUDE_DIR = ../include
SRC_DIR = ../src

# If the target is "debug", then include debugging information and build the
# marley executable with optimization turned off
ifeq ($(MAKECMDGOALS),debug)
  CXXFLAGS ?= -O0 -g -std=c++14
else
# Otherwise, use full optimization and do not include debugging info
  CXXFLAGS ?= -O3 -std=c++14
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

SHARED_LIB_NAME := MARLEY
SHARED_LIB := lib$(SHARED_LIB_NAME).$(SHARED_LIB_SUFFIX)
ROOT_SHARED_LIB_NAME := MARLEY_ROOT
ROOT_SHARED_LIB := lib$(ROOT_SHARED_LIB_NAME).$(SHARED_LIB_SUFFIX)

TEST_EXECUTABLE = martest
TEST_OBJECTS = $(notdir $(patsubst %.cc,%.o,$(wildcard $(SRC_DIR)/tests/*.cc)))

all: marley
debug: marley
test: $(TEST_EXECUTABLE)

# Skip lots of initialization if all we want is "make clean/uninstall"
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),uninstall)

  # Use g++ as the default compiler
  CXX ?= g++
  CXXFLAGS += -I$(INCLUDE_DIR) -Wall -Wextra -Wpedantic
  CXXFLAGS += -Wno-error=unused-parameter -Wcast-align

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
        CXXFLAGS += -Wno-shadow
      endif
      # Linking to ROOT libraries can be problematic on distributions (e.g.,
      # Ubuntu) that set the g++ flag -Wl,--as-needed by default (see
      # http://www.bnikolic.co.uk/blog/gnu-ld-as-needed.html for details), so
      # disable this behavior on Linux.
      ifneq ($(UNAME_S),Darwin)
        CXXFLAGS += -Wl,--no-as-needed
      endif
    endif
  endif
  
  OBJECTS := $(notdir $(patsubst %.cc,%.o,$(wildcard $(SRC_DIR)/*.cc)))
  OBJECTS := $(filter-out marley.o marley_root.o, $(OBJECTS))
  OBJECTS := $(filter-out marsum.o RootJSONConfig.o, $(OBJECTS))
  OBJECTS := $(filter-out RootOutputFile.o RootEventFileReader.o, $(OBJECTS))

  # Get information about the GNU Scientific Library installation
  # (required as of MARLEY v1.1.0)
  GSLCONFIG := $(shell which gsl-config)
  ifeq (, $(shell which gsl-config))
    define err_message
Could not find a working GNU Scientific Library installation.
Please ensure that the gsl-config executable is on your PATH and try again.
    endef
    $(error "$(err_message)")
  endif

  # Add the appropriate compilation flags to use the GNU Scientific Library
  GSL_CXXFLAGS := $(shell $(GSLCONFIG) --cflags)
  GSL_LDFLAGS := $(shell $(GSLCONFIG) --libs)

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
      CXXFLAGS += -std=c++14
    else
      ROOT_VERSION := $(shell $(ROOTCONFIG) --version)
      $(info Found ROOT version $(ROOT_VERSION) in $(ROOT))
      $(info MARLEY will be built with ROOT support.)
      override CXXFLAGS += -DUSE_ROOT
      MAYBE_MARSUM = marsum mroot
      ROOT_CXXFLAGS := $(shell $(ROOTCONFIG) --cflags)

      # If ROOT was built with C++17 support, switch to building
      # MARLEY with -std=c++17 as well. Otherwise we're likely
      # to hit problems with std::string_view in ROOT headers, etc.
      # Later -std options take precedence, so just tack on
      # a new one to keep things simple.
      ifneq (, $(findstring c++1z, $(ROOT_CXXFLAGS)))
        override CXXFLAGS += -std=c++1z
      endif
      ifneq (, $(findstring c++17, $(ROOT_CXXFLAGS)))
        override CXXFLAGS += -std=c++17
      endif
      # Do the same thing for C++2a
      ifneq (, $(findstring c++2a, $(ROOT_CXXFLAGS)))
        override CXXFLAGS += -std=c++2a
      endif

      ROOT_LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
      ROOT_LIBDIR := $(shell $(ROOTCONFIG) --libdir)
      ROOT_LDFLAGS += -L$(ROOT_LIBDIR) -lCore -lRIO -lHist -lTree
      ifeq ($(UNAME_S),Linux)
        ROOT_LDFLAGS += -rdynamic
      endif
      ROOT_DICT_INCLUDES := -I$(INCLUDE_DIR) marley/Particle.hh \
        marley/Event.hh marley/marley_linkdef.hh
      USE_ROOT = yes
      ROOT_OBJ_DICT = marley_root_dict.o

      # If we're linking against ROOT, then we need to move the
      # OutputFile.o object from the MARLEY library to the MARLEY_ROOT
      # library
      OBJECTS := $(filter-out OutputFile.o, $(OBJECTS))

      ROOT_SHARED_LIB_LDFLAGS := -l$(ROOT_SHARED_LIB_NAME)
      ROOT_SHARED_LIB_OBJECTS = marley_root.o RootJSONConfig.o
      ROOT_SHARED_LIB_OBJECTS += OutputFile.o RootOutputFile.o
      ROOT_SHARED_LIB_OBJECTS += RootEventFileReader.o $(ROOT_OBJ_DICT)
  
$(ROOT_OBJ_DICT):
	$(RM) marley_root_dict*.*
	$(ROOTCLING) -f marley_root_dict.cc -c $(ROOT_DICT_INCLUDES)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) \
	-I$(INCLUDE_DIR) -fPIC -o $(ROOT_OBJ_DICT) -c marley_root_dict.cc
	$(RM) marley_root_dict.cc
  
$(ROOT_SHARED_LIB): $(SHARED_LIB) $(ROOT_SHARED_LIB_OBJECTS)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) -fPIC -shared \
	  -o $(ROOT_SHARED_LIB) -L. -l$(SHARED_LIB_NAME) \
	  $(ROOT_SHARED_LIB_OBJECTS) $(ROOT_LDFLAGS) $(GSL_LDFLAGS)

    endif
  else
    $(info Ignoring any ROOT installations that may be present.)
    $(info MARLEY will be built without ROOT support.)
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

MARLEY_LIBS := $(SHARED_LIB)
ifeq ($(USE_ROOT),yes)
  MARLEY_LIBS += $(ROOT_SHARED_LIB)

  # If we're building the tests, then link to some extra ROOT libraries
  ifneq (,$(findstring test,$(MAKECMDGOALS)))
    ROOT_LDFLAGS += -lGraf -lGpad
  endif
endif

# Causes GNU make to auto-delete the object files when the build is complete
.INTERMEDIATE: $(OBJECTS) $(TEST_OBJECTS) $(ROOT_SHARED_LIB_OBJECTS) marley.o marsum.o

%.o: $(SRC_DIR)/%.cc
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) \
	-I$(INCLUDE_DIR) -fPIC -o $@ -c $^

%.o: $(SRC_DIR)/tests/%.cc
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) \
	-I$(INCLUDE_DIR) -fPIC -o $@ -c $^

$(SHARED_LIB): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(GSL_CXXFLAGS) $(GSL_LDFLAGS) \
	-fPIC -shared -o $@ $^

marsum: $(MARLEY_LIBS) marsum.o
	$(CXX) $(CXXFLAGS) $(GSL_CXXFLAGS) -o $@ -L. \
	  -l$(SHARED_LIB_NAME) $(ROOT_SHARED_LIB_LDFLAGS) $(ROOT_LDFLAGS) \
	  $(GSL_LDFLAGS) -Wl,-rpath -Wl,$(libdir):$(shell pwd) marsum.o

mroot: $(MARLEY_LIBS)
	cp $(SRC_DIR)/scripts/mroot .

marley: $(MARLEY_LIBS) marley.o $(MAYBE_MARSUM)
	$(CXX) $(CXXFLAGS) $(GSL_CXXFLAGS) -o $@ -L. \
	  -l$(SHARED_LIB_NAME) $(ROOT_SHARED_LIB_LDFLAGS) $(ROOT_LDFLAGS) \
	  $(GSL_LDFLAGS) -Wl,-rpath -Wl,$(libdir):$(shell pwd) marley.o

$(TEST_EXECUTABLE): $(TEST_OBJECTS) $(MARLEY_LIBS)
	$(CXX) $(CXXFLAGS) $(GSL_CXXFLAGS) -o $@ -L. \
	  -l$(SHARED_LIB_NAME) $(ROOT_SHARED_LIB_LDFLAGS) $(ROOT_LDFLAGS) \
	  $(GSL_LDFLAGS) -o $@ $(TEST_OBJECTS)

.PHONY: clean install uninstall

clean:
	$(RM) *.$(SHARED_LIB_SUFFIX) *.o marley_root_dict*.* marley
	$(RM) -rf *.dSYM marsum mroot $(TEST_EXECUTABLE)

install: marley
	mkdir -p $(DESTDIR)$(bindir)
	mkdir -p $(DESTDIR)$(libdir)
	mkdir -p $(DESTDIR)$(incdir)/marley
	mkdir -p $(DESTDIR)$(datadir)/marley
	cp marley $(MAYBE_MARSUM) $(DESTDIR)$(bindir)
	cp $(SHARED_LIB) $(DESTDIR)$(libdir)
	cp $(ROOT_SHARED_LIB) $(DESTDIR)$(libdir) 2> /dev/null || true
	cp marley_root_dict_rdict.pcm $(DESTDIR)$(libdir) 2> /dev/null || true
	cp -r ../react $(DESTDIR)$(datadir)/marley
	cp -r ../structure $(DESTDIR)$(datadir)/marley
	cp -r ../examples $(DESTDIR)$(datadir)/marley
	cp -r ../include/marley $(DESTDIR)$(incdir)
	ldconfig

uninstall:
	$(RM) $(DESTDIR)$(bindir)/marley
	$(RM) $(DESTDIR)$(bindir)/marsum
	$(RM) $(DESTDIR)$(bindir)/mroot
	$(RM) $(DESTDIR)$(libdir)/$(SHARED_LIB)
	$(RM) $(DESTDIR)$(libdir)/$(ROOT_SHARED_LIB)
	$(RM) $(DESTDIR)$(libdir)/marley_root_dict_rdict.pcm
	$(RM) -r $(DESTDIR)$(datadir)/marley
	$(RM) -r $(DESTDIR)$(incdir)/marley
	ldconfig
