#Choice of the compiler
CC=g++

#Compiling options
CFLAGS=-c

all: Mka3d

Mka3d: main.o solide.o geometry.o 
	$(CC) -o Mka3d main.o solide.o geometry.o 

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

solide.o: solide.cpp
	$(CC) $(CFLAGS) solide.cpp

geometry.o: geometry.cpp
	$(CC) $(CFLAGS) geometry.cpp

clean:
	rm *.o Mka3d mesh

mesh: mailleur.o geometry.o
	$(CC) -o mesh mailleur.o geometry.o

mailleur.o: mailleur.cpp
	$(CC) $(CFLAGS) mailleur.cpp
