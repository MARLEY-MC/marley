CXX=g++

OBJ = parse.o

%.o: %.c
	$(CXX) -c -o $@

parse: $(OBJ)
	$(CXX) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o parse
