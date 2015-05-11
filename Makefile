CXX=g++
CXXFLAGS=-std=c++11 -I. -Wall -Wextra -Werror

OBJ = react.o marley_utils.o TMarleyReaction.o #TMarleyGamma.o TMarleyLevel.o TMarleyDecayScheme.o parse.o

all: react

%.o: %.c
	$(CXX) -c -o $@

react: $(OBJ)
	$(CXX) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o react
