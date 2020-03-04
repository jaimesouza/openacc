# Your compiler
CC = pgcc

# flags
CFLAGS = -fast -ta=turing -Minfo=all
#CFLAGS =

default: all

all: serial gpu

serial:
	$(CC) $(CFLAGS) -o serial-jacobi-laplace.o jacobi-laplace-serial.c

gpu:
	$(CC) $(CFLAGS) -o openacc-jacobi-laplace.o jacobi-laplace-openacc.c

clean:
	rm -f *.o
.PHONY: clean
