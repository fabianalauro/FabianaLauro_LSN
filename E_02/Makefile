CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: E_02_1.exe E_02_2.exe 

E_02_1.exe : E_02_1.o random.o
	$(CC) E_02_1.o random.o -o E_02_1.exe
E_02_1.o : E_02_1.cpp random.h
	$(CC) -c E_02_1.cpp -o E_02_1.o $(CFLAGS)
	
E_02_2.exe : E_02_2.o random.o Position.o
	$(CC) E_02_2.o random.o Position.o -o E_02_2.exe
E_02_2.o : E_02_2.cpp random.h Position.hpp
	$(CC) -c E_02_2.cpp -o E_02_2.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	
Position.o : Position.cpp Position.hpp
	$(CC) -c Position.cpp -o Position.o $(CFLAGS)
	
Fun_01.o : Fun_01.cpp Fun_01.h
	$(CC) -c Fun_01.cpp -o Fun_01.o $(CFLAGS)

clean :
	rm *.o *.exe *.dat *.out
	
	
