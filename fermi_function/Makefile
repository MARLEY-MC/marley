CXX=g++
CXXFLAGS=-std=c++11 -I. -Wall -Wextra -Werror

OBJ =  marley_utils.o fermi_function.o

%.o: %.c
	$(CXX) -c -o $@

fermi: $(OBJ)
	$(CXX) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o fermi
