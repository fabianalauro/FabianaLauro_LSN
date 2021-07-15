CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: E_09.exe 

E_09.exe : E_09.o random.o Objects.o 
	$(CC) E_09.o random.o Objects.o -o E_09.exe
	
E_09.o : E_09.cpp random.h Objects.hpp
	$(CC) -c E_09.cpp -o E_09.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	
Objects.o : Objects.cpp Objects.hpp random.h 
	$(CC) -c Objects.cpp -o Objects.o $(CFLAGS)
	
clean :
	rm *.o *.exe *.dat *.out
	
