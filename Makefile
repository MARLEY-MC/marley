CXX=g++
CXXFLAGS=-std=c++11 -I. -Wall -Wextra -Wpedantic -Werror #-O3
USE_ROOT=yes

OBJ = marley_utils.o TMarleyParticle.o TMarleyEvent.o TMarleyEvaporationThreshold.o
OBJ += TMarleyReaction.o TMarleyGamma.o TMarleyLevel.o TMarleyDecayScheme.o
OBJ += TMarleyMassTable.o

ifdef USE_ROOT
# Adding the g++ compiler option -DUSE_ROOT to the CXXFLAGS
# variable allows you to use conditional compilation
# via the preprocessor directive #ifdef USE_ROOT.
# Currently none of the core MARLEY classes use such
# preprocessor directives, but the example executable
# react does.
#
# The root_dict.o object file should be added
# to the list of prerequisites for any executable
# that uses TTrees containing TMarleyEvents
# or TMarleyParticles
CXXFLAGS += `root-config --cflags` -DUSE_ROOT
OBJ_DICT = root_dict.o
LDFLAGS=`root-config --libs`
endif

all: parse react validate

%.o: %.c
	$(CXX) -c -o $@

parse: $(OBJ) parse.o
	$(CXX) -o $@ $^

react: $(OBJ) $(OBJ_DICT) react.o
	$(CXX) -o $@ $^ $(LDFLAGS)

validate: $(OBJ) $(OBJ_DICT) validate.o
	$(CXX) -o $@ $^ $(LDFLAGS)

# Add more header files to the prerequisites for
# root_dict.o if you would like to store other
# MARLEY classes in ROOT trees. All such classes
# currently use a single automatically-generated
# dictionary source file root_dict.cc
# 
# The commands invoked to create root_dict.o
# do the following things:
# 1. Remove old dictionary files
# 2. Create new dictionary files, enabling ROOT
#    i/o by adding the '+' suffix to each prerequisite
#    header file (.hh extension)
# 3. Compile the dictionary source file
root_dict.o: TMarleyParticle.hh TMarleyEvent.hh
	rm -f root_dict.cc root_dict.h
	rootcint root_dict.cc -c $(subst .hh,.hh+,$^)
	$(CXX) $(CXXFLAGS) -c -o root_dict.o root_dict.cc

.PHONY: clean

clean:
	rm -f *.o parse react validate root_dict.cc root_dict.h
