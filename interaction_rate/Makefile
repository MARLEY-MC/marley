CXX=g++
CXXFLAGS=-std=c++11 -I. -Wall -Wextra -Werror

rate_spectra: rate_spectra.cc
	$(CXX) -o rate_spectra rate_spectra.cc

.PHONY: clean

clean:
	rm -f *.o rate_spectra
