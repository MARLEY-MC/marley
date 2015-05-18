CXX=g++
CXXFLAGS=-std=c++11 -I.  -Wall -Wextra `root-config --cflags`
OBJ = threshold.o

all: threshold

%.o: %.c
	$(CXX) -c -o $@ `root-config --libs`

threshold: $(OBJ)
	$(CXX) -o $@ $^ `root-config --libs`


.PHONY: clean

clean:
	rm -f *.o threshold out.txt energy_hist.root
