CXX=g++
CXXFLAGS=-std=c++11 -O3 #-Wall -Wextra -Wpedantic -Werror -Wno-error=unused-parameter

OBJ = optmod.o /home/sjg/Desktop/marley/TMarleyMassTable.o

all: opt

%.o: %.c
	$(CXX) -c -o $@

opt: $(OBJ)
	$(CXX) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o opt
