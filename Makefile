CXX=g++
CXXFLAGS=-std=c++11 -I. -Wall -Wextra -Werror

OBJ = TMarleyGamma.o TMarleyLevel.o TMarleyDecayScheme.o marley_utils.o parse.o

%.o: %.c
	$(CXX) -c -o $@

parse: $(OBJ)
	$(CXX) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o parse
