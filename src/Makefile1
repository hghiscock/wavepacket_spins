CC = g++

all: test

test: procedures.o hamiltonians.o param.o results.o sc_tensor.o calcphis.o
	$(CC) -Wall procedures.o param.o hamiltonians.o results.o sc_tensor.o calcphis.o -o test

calcphis.o: calcphis.cpp
	$(CC) -c calcphis.cpp

results.o: results.cpp
	$(CC) -c results.cpp

sc_tensor.o: sc_tensor.cpp
	$(CC) -c sc_tensor.cpp

hamiltonians.o: hamiltonians.cpp
	$(CC) -c hamiltonians.cpp

param.o: param.cpp
	$(CC) -c param.cpp 

procedures.o: procedures.cpp
	$(CC) -c procedures.cpp 

clean:
	rm *.o test
