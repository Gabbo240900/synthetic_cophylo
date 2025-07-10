CC = gcc
CCFLAGS = -O3 -lm

cophylo: sources/main_cophylo.o
	$(CC) -o bin/cophylo.out sources/main_cophylo.o $(CCFLAGS)

main_cophylo.o: sources/main_cophylo.c 
	$(CC) -c -g sources/main_cophylo.c $(CCFLAGS)

clean:
	rm sources/*.o
