INCLUDE_DIR = ..
SRC_DIR = ..

ifndef MARLEY
  $(error The MARLEY environment variable is not set. Please source the \
     setup_marley.sh script and try again.)
else
  MARCONF := $(shell command -v marley-config 2> /dev/null)
  ifndef MARCONF
    $(error The marley-config script is not present on the system PATH. \
        Please ensure that MARLEY has been fully built and try again.)
  endif
endif

# If the target is "debug", then include debugging information and build the
# example programs with optimization turned off
ifeq ($(MAKECMDGOALS),debug)
  CXXFLAGS ?= -O0 -g
else
# Otherwise, use full optimization and do not include debugging info
  CXXFLAGS ?= -O3
endif

# Use GCC as the default compiler. Enable a few useful warnings.
CXX = g++
CXXFLAGS += -Wall -Wextra -Wpedantic -Wcast-align

all: marprint mardumpxs
debug: all

# Use the marley-config script to get the MARLEY compiler flags and
# library names.
override CXXFLAGS += $(shell marley-config --cflags)
MARLEY_LIBS = $(shell marley-config --libs)

OBJECTS := $(notdir $(patsubst %.cc,%.o,$(wildcard $(SRC_DIR)/*.cc)))

.INTERMEDIATE: $(OBJECTS)

%.o: $(SRC_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -o $@ -c $^

marprint: marprint.o
	$(CXX) $(CXXFLAGS) -o $@ $(MARLEY_LIBS) marprint.o

mardumpxs: mardumpxs.o
	$(CXX) $(CXXFLAGS) -o $@ $(MARLEY_LIBS) mardumpxs.o

mardecay: mardecay.o
	$(CXX) $(CXXFLAGS) -o $@ $(MARLEY_LIBS) mardecay.o

.PHONY: clean

clean:
	$(RM) *.o marprint mardumpxs mardecay
