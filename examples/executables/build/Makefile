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

all: marprint martxs mardxs marjoin
debug: all

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

martxs: martxs.o
	$(CXX) $(CXXFLAGS) -o $@ $(MARLEY_LIBS) $(ROOT_LIBS) martxs.o

mardxs: mardxs.o
	$(CXX) $(CXXFLAGS) -o $@ $(MARLEY_LIBS) $(ROOT_LIBS) mardxs.o

marjoin: marjoin.o
	$(CXX) $(CXXFLAGS) -o $@ $(MARLEY_LIBS) $(ROOT_LIBS) marjoin.o

.PHONY: clean

clean:
	$(RM) *.o marprint martxs mardxs marjoin
