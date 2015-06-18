CXX=g++
CXXFLAGS=-std=c++11 -I. -Wall -Wextra -Werror
USE_ROOT=yes

OBJ = marley_utils.o TMarleyParticle.o TMarleyEvent.o TMarleyEvaporationThreshold.o TMarleyReaction.o TMarleyGamma.o TMarleyLevel.o TMarleyDecayScheme.o TMarleyMassTable.o react.o #parse.o

ifdef USE_ROOT
# Adding the g++ compiler option -DUSE_ROOT to the CXXFLAGS
# variable allows you to use conditional compilation
# via the preprocessor directive #ifdef USE_ROOT.
# Currently none of the core MARLEY classes use such
# preprocessor directives, but the example executable
# react does.
CXXFLAGS += `root-config --cflags` -DUSE_ROOT
OBJ += root_dict.o
LDFLAGS=`root-config --libs`
endif

all: react

%.o: %.c
	$(CXX) -c -o $@

react: $(OBJ)
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

#parse: $(OBJ)
#	$(CXX) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o react root_dict.cc root_dict.h
