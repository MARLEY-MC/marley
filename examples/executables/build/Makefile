INCLUDE_DIR = ..
SRC_DIR = ..

ifndef MARLEY
  $(error The MARLEY environment variable is not set. Please source the \
     setup_marley.sh script and try again.)
endif

MARLEY_INC = $(MARLEY)/include
MARLEY_LIB = $(MARLEY)/build

# If the target is "debug", then include debugging information and build the
# example programs with optimization turned off
ifeq ($(MAKECMDGOALS),debug)
  CXXFLAGS ?= -O0 -g -std=c++14
else
# Otherwise, use full optimization and do not include debugging info
  CXXFLAGS ?= -O3 -std=c++14
endif

CXX = g++
CXXFLAGS += -Wall -Wextra -Wpedantic -Wcast-align

all: marprint mardumpxs
debug: all

# Always link against the MARLEY shared library. If the MARLEY_ROOT shared
# library is present in the build/ folder, then link against it. Also enable
# ROOT-dependent features by defining the USE_ROOT preprocessor macro.
MARLEY_LIBS := -L$(MARLEY_LIB) -lMARLEY
ifneq (,$(wildcard $(MARLEY)/build/libMARLEY_ROOT*))
  MARLEY_LIBS += -lMARLEY_ROOT
  override CXXFLAGS += -DUSE_ROOT $(shell root-config --cflags)
  ROOT_LIBS = $(shell root-config --libs)
endif

OBJECTS := $(notdir $(patsubst %.cc,%.o,$(wildcard $(SRC_DIR)/*.cc)))

.INTERMEDIATE: $(OBJECTS)

%.o: $(SRC_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -I$(MARLEY_INC) -o $@ -c $^

marprint: marprint.o
	$(CXX) $(CXXFLAGS) -o $@ $(MARLEY_LIBS) $(ROOT_LIBS) marprint.o

mardumpxs: mardumpxs.o
	$(CXX) $(CXXFLAGS) -o $@ $(MARLEY_LIBS) $(ROOT_LIBS) mardumpxs.o

.PHONY: clean

clean:
	$(RM) *.o marprint mardumpxs
