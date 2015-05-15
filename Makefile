CXX=g++
CXXFLAGS=-std=c++11 -I. -Wall -Wextra -Werror

OBJ = marley_utils.o TMarleyReaction.o TMarleyGamma.o TMarleyLevel.o TMarleyDecayScheme.o parse.o #react.o 

all: parse

%.o: %.c
	$(CXX) -c -o $@

#react: $(OBJ)
#	$(CXX) -o $@ $^

parse: $(OBJ)
	$(CXX) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o react
