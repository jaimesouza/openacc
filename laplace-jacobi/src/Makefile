# Your compiler
CC = pgcc

# flags
CFLAGS = -fast -acc -ta=nvidia,time -Minfo=all

default: all

all: serial gpu

serial:
	$(CC) $(CFLAGS) -o serial-jacobi-laplace.o jacobi-laplace-serial.c

gpu:
	$(CC) $(CFLAGS) -o openacc-jacobi-laplace.o jacobi-laplace-openacc.c

clean:
	rm -f *.o
.PHONY: clean
