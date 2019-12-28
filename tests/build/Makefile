TEST_EXECUTABLE = test_marley

INCLUDE_DIR = $(MARLEY)/include
LIB_DIR = $(MARLEY)/lib
TESTS_DIR = $(MARLEY)/tests

CXXFLAGS ?= -O3 -std=c++14

all: $(TEST_EXECUTABLE)

# TODO: reduce code duplication between this Makefile and the main one
ifneq ($(MAKECMDGOALS),clean)

  ifndef MARLEY
    $(error The MARLEY environment variable is not set)
  endif

  # Use g++ as the default compiler
  CXX ?= g++
  CXXFLAGS += -I$(INCLUDE_DIR) -Wall -Wextra -Wpedantic
  CXXFLAGS += -Wno-error=unused-parameter -Wcast-align

  OBJECTS := $(notdir $(patsubst %.cc,%.o,$(wildcard $(TESTS_DIR)/src/*.cc)))

  # Get information about the GNU Scientific Library installation
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
      ROOT_LDFLAGS += -L$(ROOT_LIBDIR) -lCore -lHist -lTree -lGraf -lGpad
      ifeq ($(UNAME_S),Linux)
        ROOT_LDFLAGS += -rdynamic
      endif
      USE_ROOT = yes
    endif
  else
    $(info Ignoring any ROOT installations that may be present.)
    $(info MARLEY tests will be built without ROOT support.)
  endif
endif

MARLEY_LIBS := -lMARLEY
ifeq ($(USE_ROOT),yes)
  MARLEY_LIBS += -lMARLEY_ROOT
endif

# Causes GNU make to auto-delete the object files when the build is complete
.INTERMEDIATE: $(OBJECTS)

%.o: $(TESTS_DIR)/src/%.cc
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) \
	-I$(INCLUDE_DIR) -fPIC -o $@ -c $^

$(TEST_EXECUTABLE): $(OBJECTS)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) $(GSL_LDFLAGS) \
	$(ROOT_LDFLAGS) $(MARLEY_LIBS) -fPIC -o $@ $^

.PHONY: clean

clean:
	$(RM) *.o $(TEST_EXECUTABLE) *.pdf
