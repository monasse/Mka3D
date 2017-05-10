#Choice of the compiler
CC=g++

#Compiling options
CFLAGS=-c

all: Mka3d

Mka3d: main.o solide.o geometry.o forces_ext.o
	$(CC) -o Mka3d main.o solide.o geometry.o forces_ext.o

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

solide.o: solide.cpp
	$(CC) $(CFLAGS) solide.cpp

geometry.o: geometry.cpp
	$(CC) $(CFLAGS) geometry.cpp

forces_ext.o: forces_ext.cpp
	$(CC) $(CFLAGS) forces_ext.cpp

vitesse.o: vitesse.cpp
	$(CC) $(CFLAGS) vitesse.cpp

clean:
	rm *.o Mka3d mesh

mesh: mailleur.o geometry.o vitesse.o
	$(CC) -o mesh mailleur.o geometry.o vitesse.o

mailleur.o: mailleur.cpp
	$(CC) $(CFLAGS) mailleur.cpp
