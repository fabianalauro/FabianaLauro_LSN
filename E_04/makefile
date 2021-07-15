CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: MolDyn_NVE.exe

MolDyn_NVE.exe : MolDyn_NVE.cpp MolDyn_NVE.h BlockingAve.o
	$(CC) MolDyn_NVE.cpp BlockingAve.o -o MolDyn_NVE.exe
	
BlockingAve.o : BlockingAve.cpp BlockingAve.hpp
	$(CC) -c BlockingAve.cpp -o BlockingAve.o $(CFLAGS)

clean :
	rm *.o rm *.exe
	
	
