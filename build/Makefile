INCLUDE_DIR = ../include
SRC_DIR = ../src

# Set C++ compiler flags common to all build types
CXXFLAGS = -std=c++14 -I$(INCLUDE_DIR) -Wall -Wextra -Wpedantic
CXXFLAGS += -Werror -Wno-error=unused-parameter -Wcast-align

# This Makefile provides three targets with different debugging
# and optimization compiler flags. The default, "prof", enables
# optimization but includes debugging information. The "debug"
# target turns optimization off and, unsurprisingly, also includes
# debugging information. The third, "opt", enables optimization
# and omits debugging information.
prof: CXXFLAGS += -O3 -g -fno-omit-frame-pointer -DNDEBUG
prof: all

debug: CXXFLAGS += -O0 -g
debug: all

opt: CXXFLAGS += -O3 -DNDEBUG
opt: all

all: marley $(SHARED_LIB) $(ROOT_SHARED_LIB)

# Skip lots of initialization if all we want is "make clean/install/uninstall"
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),install)
ifneq ($(MAKECMDGOALS),uninstall)
 
  # Use g++ as the default compiler
  CXX = g++

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
    ifneq (,$(or $(findstring GCC,$(CXXVERSION)), \
      $(findstring g++,$(CXXVERSION))))

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
      # disable this behavior
      CXXFLAGS += -Wl,--no-as-needed
    endif
  endif
  
  OBJECTS := $(notdir $(patsubst %.cc,%.o,$(wildcard $(SRC_DIR)/*.cc)))
  OBJECTS := $(filter-out marley.o RootConfigurationFile.o, $(OBJECTS))
  
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
    USE_ROOT=""
  else
    ROOT_VERSION := $(shell $(ROOTCONFIG) --version)
    $(info Found ROOT version $(ROOT_VERSION) in $(ROOT))
    $(info MARLEY will be built with ROOT support.)
    CXXFLAGS += -DUSE_ROOT
    ROOT_CXXFLAGS := $(shell $(ROOTCONFIG) --cflags)
    ROOT_LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
    ROOT_LIBDIR := $(shell $(ROOTCONFIG) --libdir)
    ROOT_LDFLAGS := $(ROOT_LDFLAGS) -L$(ROOT_LIBDIR) -lCore -lRIO \
      -lHist -lTree -rdynamic
    ROOT_DICT_INCLUDES := -I$(INCLUDE_DIR) marley/Particle.hh \
      marley/Event.hh marley/marley_linkdef.hh
    USE_ROOT = yes
    ROOT_OBJ_DICT = marley_root_dict.o
    ROOT_SHARED_LIB_NAME := MARLEY_ROOT
    ROOT_SHARED_LIB_LDFLAGS := -l$(ROOT_SHARED_LIB_NAME)
    ROOT_SHARED_LIB = lib$(ROOT_SHARED_LIB_NAME).$(SHARED_LIB_SUFFIX)
    ROOT_SHARED_LIB_OBJECTS = RootConfigurationFile.o $(ROOT_OBJ_DICT)
  
$(ROOT_OBJ_DICT):
	rm -f marley_root_dict*.*
	$(ROOTCLING) -f marley_root_dict.cc -c $(ROOT_DICT_INCLUDES)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) -fPIC -o $(ROOT_OBJ_DICT) -c \
	  marley_root_dict.cc
  
$(ROOT_SHARED_LIB): $(SHARED_LIB) $(ROOT_SHARED_LIB_OBJECTS)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) -fPIC -shared \
	  -o $(ROOT_SHARED_LIB) -L. -l$(SHARED_LIB_NAME) \
	  $(ROOT_SHARED_LIB_OBJECTS) $(ROOT_LD_FLAGS)

  endif
endif
endif
endif

# Set up default variables for the install/uninstall targets and for setting
# the executable rpath
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
datadir = $(prefix)/share
libdir = $(exec_prefix)/lib
incdir = $(prefix)/include

%.o: $(SRC_DIR)/%.cc
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) -fPIC -o $@ -c $^

$(SHARED_LIB): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -fPIC -shared -o \
	  lib$(SHARED_LIB_NAME).$(SHARED_LIB_SUFFIX) $^

marley: $(SHARED_LIB) $(ROOT_SHARED_LIB) marley.o
	$(CXX) $(CXXFLAGS) -o $@ -L. \
	  -l$(SHARED_LIB_NAME) $(ROOT_SHARED_LIB_LDFLAGS) $(ROOT_LDFLAGS) \
	  -Wl,-rpath -Wl,$(libdir):$(shell pwd) marley.o

.PHONY: clean install uninstall

clean:
	rm -f *.so *.o marley_root_dict*.* marley

install:
	mkdir -p $(DESTDIR)$(bindir)
	mkdir -p $(DESTDIR)$(libdir)
	mkdir -p $(DESTDIR)$(incdir)/marley
	mkdir -p $(DESTDIR)$(datadir)/marley
	cp marley $(DESTDIR)$(bindir)
	cp libMARLEY* $(DESTDIR)$(libdir)
	cp -r ../react $(DESTDIR)$(datadir)/marley
	cp -r ../examples $(DESTDIR)$(datadir)/marley
	cp -r ../structure $(DESTDIR)$(datadir)/marley
	cp -r ../include/marley $(DESTDIR)$(incdir)

uninstall:
	rm -f $(DESTDIR)$(bindir)/marley
	rm -f $(DESTDIR)$(libdir)/libMARLEY*
	rm -rf $(DESTDIR)$(datadir)/marley
	rm -rf $(DESTDIR)$(incdir)/marley
