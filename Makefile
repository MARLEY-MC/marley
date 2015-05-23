CXX=g++
CXXFLAGS=-std=c++11 -I. -Wall -Wextra #-Werror

OBJ = marley_utils.o TMarleyParticle.o TMarleyEvent.o TMarleyReaction.o TMarleyGamma.o TMarleyLevel.o TMarleyDecayScheme.o react.o #parse.o

all: react

%.o: %.c
	$(CXX) -c -o $@

react: $(OBJ)
	$(CXX) -o $@ $^

#parse: $(OBJ)
#	$(CXX) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o react
