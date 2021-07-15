CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: E_05.exe

E_05.exe : E_05.o random.o Metropolis.o Probability.o BlockingAve.o
	$(CC) E_05.o random.o Metropolis.o Probability.o BlockingAve.o -o E_05.exe
E_05.o : E_05.cpp random.h Metropolis.hpp Probability.hpp BlockingAve.hpp
	$(CC) -c E_05.cpp -o E_05.o $(CFLAGS)


random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	
Metropolis.o : Metropolis.cpp Metropolis.hpp random.h Probability.hpp
	$(CC) -c Metropolis.cpp -o Metropolis.o $(CFLAGS)
	
Probability.o : Probability.cpp Probability.hpp
	$(CC) -c Probability.cpp -o Probability.o $(CFLAGS)
	
BlockingAve.o : BlockingAve.cpp BlockingAve.hpp
	$(CC) -c BlockingAve.cpp -o BlockingAve.o $(CFLAGS)
	
clean :
	rm *.o *.exe *.dat *.out
	
