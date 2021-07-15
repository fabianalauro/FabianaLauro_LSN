CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: E_10.1.exe E_10.2.exe

E_10.1.exe : E_10.1.o random.o Objects.o 
	mpicxx E_10.1.o random.o Objects.o -o E_10.1.exe
	
E_10.1.o : E_10.1.cpp random.h Objects.hpp
	mpicxx -c E_10.1.cpp -o E_10.1.o $(CFLAGS)

E_10.2.exe : E_10.2.o random.o Objects.o 
	mpicxx E_10.2.o random.o Objects.o -o E_10.2.exe

E_10.2.o : E_10.2.cpp random.h Objects.hpp
	mpicxx -c E_10.2.cpp -o E_10.2.o $(CFLAGS)

random.o : random.cpp random.h
	mpicxx -c random.cpp -o random.o $(CFLAGS)
	
Objects.o : Objects.cpp Objects.hpp random.h 
	mpicxx -c Objects.cpp -o Objects.o $(CFLAGS)
	
clean :
	rm *.o *.exe *.dat *.out
	

#To run E_10.2.exe: mpiexec --hostfile hostfile -np 4 E_10.2.exe


	
