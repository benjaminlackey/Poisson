#Specify compiler
CC=gcc

#Compiling options
CFLAGS=-g -lm -lgsl -lgslcblas -Wall -pedantic -ansi -O2 -W

OBJ=main.o print.o coefficients.o matrixoperators.o integratesphere.o poisson.o

#Main compiling instructions 
poisson: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o poisson

main.o: poisson.h main.c
	$(CC) -c $(CFLAGS) main.c

print.o: poisson.h print.c
	$(CC) -c $(CFLAGS) print.c

coefficients.o: poisson.h coefficients.c
	$(CC) -c $(CFLAGS) coefficients.c

matrixoperators.o: poisson.h matrixoperators.c
	$(CC) -c $(CFLAGS) matrixoperators.c

integratesphere.o: poisson.h integratesphere.c
	$(CC) -c $(CFLAGS) integratesphere.c

poisson.o: poisson.h poisson.c
	$(CC) -c $(CFLAGS) poisson.c
