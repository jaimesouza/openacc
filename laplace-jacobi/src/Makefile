# Your compiler
CC = pgcc

# flags
CFLAGS = -fast -ta=turing -Minfo=all
#CFLAGS =

default: all

all:
	$(CC) $(CFLAGS) -o jacobi-laplace-serial.o jacobi-laplace-serial.c

clean:
	rm -f *.o
.PHONY: clean
