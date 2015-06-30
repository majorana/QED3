CC = gcc -std=c99 -pedantic -Wall -O3

all: qed

fields.o: fields.c fields.h lattice.h linalg.h complex/complex.h rand/ranlxd.h Makefile
	$(CC) -c $< -o $@

hmc.o: hmc.c hmc.h lattice.h fields.h linalg.h complex/complex.h rand/ranlxd.h rand/gauss.h Makefile
	$(CC) -c $< -o $@

lattice.o: lattice.c lattice.h
	$(CC) -c $< -o $@

linalg.o: linalg.c linalg.h lattice.h complex/complex.h Makefile
	$(CC) -c $< -o $@

ranlxd.o: rand/ranlxd.c rand/ranlxd.h Makefile
	$(CC) -c $< -o $@ -I rand

gauss.o: rand/gauss.c rand/gauss.h Makefile
	$(CC) -c $< -o $@ -I rand

fermion.o: fermion.c fermion.h lattice.h linalg.h complex/complex.h rand/ranlxd.h Makefile
	$(CC) -c $< -o $@

integrator.o: integrator.c integrator.h hmc.h Makefile
	$(CC) -c $< -o $@

qed.o: qed.c fields.h lattice.h linalg.h hmc.h complex/complex.h Makefile rand/ranlxd.h
	$(CC) -c $< -o $@

qed: fields.o qed.o hmc.o lattice.o linalg.o measurements.o ranlxd.o gauss.o fermion.o integrator.o Makefile
	$(CC) qed.o fields.o hmc.o lattice.o linalg.o measurements.o ranlxd.o gauss.o fermion.o integrator.o -o qed -lm

clean:
	rm -f *.o qed test1 test2
