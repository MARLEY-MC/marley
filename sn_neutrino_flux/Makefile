CXX=g++
CXXFLAGS=-std=c++11 -I. -Wall -Wextra -Werror

sn_neutrino_spectra: sn_neutrino_spectra.cc
	$(CXX) -o sn_neutrino_spectra sn_neutrino_spectra.cc

.PHONY: clean

clean:
	rm -f *.o sn_neutrino_spectra
