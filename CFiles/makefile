OBJS	= AuxillaryFunctions.o Graph.o MT19937.o Driver.o
SOURCE	= AuxillaryFunctions.c Graph.c MT19937.c Driver.c
HEADER	= AuxillaryFunctions.h Graph.h MT19937.h
OUT	= Main
CC	 = gcc -no-pie -pthread
FLAGS	 = -g -c -Wall -lm
LFLAGS	 = 

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

AuxillaryFunctions.o: AuxillaryFunctions.c
	$(CC) $(FLAGS) AuxillaryFunctions.c 

Graph.o: Graph.c
	$(CC) $(FLAGS) Graph.c 

MT19937.o: MT19937.c
	$(CC) $(FLAGS) MT19937.c 

Driver.o: Driver.c
	$(CC) $(FLAGS) Driver.c 


clean:
	rm -f $(OBJS) $(OUT)
