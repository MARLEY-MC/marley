CXX=g++
CXXFLAGS=-g -O3 -std=c++14 -I. -fPIC -Wall -Wextra -Wpedantic -Werror -Wno-error=unused-parameter
USE_ROOT=yes

OBJ = marley_utils.o meta_numerics.o Particle.o Event.o
OBJ += Generator.o Reaction.o NuclearReaction.o
OBJ += ElectronReaction.o Gamma.o Level.o
OBJ += DecayScheme.o MassTable.o StructureDatabase.o
OBJ += ConfigFile.o NuclearPhysics.o
OBJ += BackshiftedFermiGasModel.o SphericalOpticalModel.o
OBJ += NeutrinoSource.o Kinematics.o DecayChannel.o
OBJ += Integrator.o

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
# that uses TTrees containing Events
# or Particles
CXXFLAGS += `root-config --cflags` -DUSE_ROOT
OBJ_DICT = marley_root_dict.o
LDFLAGS=`root-config --libs`
endif

all: marley sharedlib

%.o: %.c
	$(CXX) -c -o $@

marley: $(OBJ) $(OBJ_DICT) marley.o
	$(CXX) -o $@ $^ $(LDFLAGS)

sharedlib: $(OBJ) $(OBJ_DICT)
	$(CXX) -shared -o libMARLEY.so $^ $(LDFLAGS)

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
marley_root_dict.o: Particle.hh Event.hh marley_linkdef.hh
	rm -f marley_root_dict.cc marley_root_dict.h
	rootcint marley_root_dict.cc -c $^
	$(CXX) $(CXXFLAGS) -c -o marley_root_dict.o marley_root_dict.cc

.PHONY: clean

clean:
	rm -f *.o *.so marley marley_root_dict.*
