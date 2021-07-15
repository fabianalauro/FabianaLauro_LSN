CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: E_01_1.exe E_01_2.exe E_01_3.exe 

E_01_1.exe : E_01_1.o random.o
	$(CC) E_01_1.o random.o -o E_01_1.exe
E_01_1.o : E_01_1.cpp random.h
	$(CC) -c E_01_1.cpp -o E_01_1.o $(CFLAGS)
	
E_01_2.exe : E_01_2.o random.o
	$(CC) E_01_2.o random.o -o E_01_2.exe
E_01_2.o : E_01_2.cpp random.h
	$(CC) -c E_01_2.cpp -o E_01_2.o $(CFLAGS)
	
E_01_3.exe : E_01_3.o random.o Fun_01.o
	$(CC) E_01_3.o random.o Fun_01.o -o E_01_3.exe
E_01_3.o : E_01_3.cpp random.h Fun_01.h
	$(CC) -c E_01_3.cpp -o E_01_3.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	
Fun_01.o : Fun_01.cpp Fun_01.h
	$(CC) -c Fun_01.cpp -o Fun_01.o $(CFLAGS)

clean :
	rm *.o *.exe *.dat *.out
	
	
