CXX=g++
CXXFLAGS=-std=c++11 -I. -Wall -Wextra -Werror `root-config --cflags`

OBJ = marley_utils.o TMarleyParticle.o TMarleyEvent.o TMarleyReaction.o TMarleyGamma.o TMarleyLevel.o TMarleyDecayScheme.o react.o #parse.o

all: react

%.o: %.c
	$(CXX) -c -o $@

react: $(OBJ)
	$(CXX) -o $@ $^ `root-config --libs`

#parse: $(OBJ)
#	$(CXX) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o react
