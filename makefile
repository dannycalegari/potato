CC=g++
CFLAGS=-g -Wall -O3 -fast
IFLAGS=-I/usr/X11R6/include
LFLAGS=-L/usr/X11R6/lib -lX11
all: potato

potato: potato.cc lobachevsky.cc vector.cc graphics.cc packing.cc input_output.cc dual_format.cc layout.cc
	$(CC) $(CFLAGS) $(IFLAGS) -o potato potato.cc $(LFLAGS) -lm

clean: rm potato
