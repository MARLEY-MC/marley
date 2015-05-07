CXX=g++
CXXFLAGS=-std=c++11 -I. -Wall -Wextra -Werror

OBJ = TMarleyGamma.o TMarleyLevel.o TMarleyDecayScheme.o marley_utils.o parse.o

all: parse sn_neutrino_spectra

%.o: %.c
	$(CXX) -c -o $@

parse: $(OBJ)
	$(CXX) -o $@ $^

sn_spectra: sn_neutrino_spectra.cc
	$(CXX) -o sn_neutrino_spectra sn_neutrino_spectra.cc

.PHONY: clean

clean:
	rm -f *.o parse sn_neutrino_spectra
