CC = g++

CFLAGS = -pedantic -W -Wall -Wextra -O3 -mtune=native

ALLCPP = CDForce.cpp
ALLHPP = CDForce.hpp

all: run

run : $(ALLCPP) $(ALLHPP) Makefile
	$(CC) $(CFLAGS) -o run CDForce.cpp

clean:
	rm -f *.o run 

