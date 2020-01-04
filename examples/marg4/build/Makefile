INCLUDE_DIR = ../include
SRC_DIR = ../src

ifndef MARLEY
  $(error The MARLEY environment variable is not set. Please source the \
     setup_marley.sh script and try again.)
endif

MARLEY_INC = $(MARLEY)/include
MARLEY_LIB = $(MARLEY)/build

# If the target is "debug", then include debugging information and build the
# Geant4 example with optimization turned off
ifeq ($(MAKECMDGOALS),debug)
  CXXFLAGS ?= -O0 -g -std=c++14
else
# Otherwise, use full optimization and do not include debugging info
  CXXFLAGS ?= -O3 -std=c++14
endif

CXX = g++
CXXFLAGS += -Wall -Wextra -Wpedantic -Wcast-align

all: marg4
debug: marg4

G4_CXXFLAGS := $(shell geant4-config --cflags)

# If Geant4 was built with C++17 or C++2a support, switch to building the
# example with a more up-to-date C++ standard. Later -std options take
# precedence, so just tack on a new one to keep things simple.
ifneq (, $(findstring c++1z, $(G4_CXXFLAGS)))
  override CXXFLAGS += -std=c++1z
endif
ifneq (, $(findstring c++17, $(G4_CXXFLAGS)))
  override CXXFLAGS += -std=c++17
endif
ifneq (, $(findstring c++2a, $(G4_CXXFLAGS)))
  override CXXFLAGS += -std=c++2a
endif

MARLEY_LIBS := -L$(MARLEY_LIB) -lMARLEY
ifneq (,$(wildcard $(MARLEY)/build/libMARLEY_ROOT*))
  MARLEY_LIBS += -lMARLEY_ROOT
  override CXXFLAGS += -DUSE_ROOT
endif

G4_LIBS := $(shell geant4-config --libs)

OBJECTS := $(notdir $(patsubst %.cc,%.o,$(wildcard $(SRC_DIR)/*.cc)))

.INTERMEDIATE: $(OBJECTS) marg4.o

%.o: $(SRC_DIR)/%.cc
	$(CXX) $(G4_CXXFLAGS) $(CXXFLAGS) -I$(INCLUDE_DIR) \
	  -I$(MARLEY_INC) -o $@ -c $^

marg4: $(OBJECTS)
	$(CXX) $(G4_CXXFLAGS) $(CXXFLAGS) -o $@ $(MARLEY_LIBS) \
	  $(G4_LIBS) $(OBJECTS)

.PHONY: clean

clean:
	$(RM) *.o marg4
