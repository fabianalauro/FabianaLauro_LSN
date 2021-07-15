CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: E_03.exe

E_03.exe : E_03.o random.o PlainVanilla.o
	$(CC) E_03.o random.o PlainVanilla.o -o E_03.exe
E_02_1.o : E_03.cpp random.h PlainVanilla.hpp
	$(CC) -c E_03.cpp -o E_03.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	
PlainVanilla.o : PlainVanilla.cpp PlainVanilla.hpp
	$(CC) -c PlainVanilla.cpp -o PlainVanilla.o $(CFLAGS)

clean :
	rm *.o *.exe *.dat *.out
	
	
