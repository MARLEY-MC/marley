CXX=g++
CXXFLAGS=-std=c++11

OBJ = parse.o

%.o: %.c
	$(CXX) -c -o $@

parse: $(OBJ)
	$(CXX) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o parse
