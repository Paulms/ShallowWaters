# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include template


.DEFAULT:
	-echo $@ does not exist.
all: ../bin/SWSolver
RiemannSolvers.o: ./RiemannSolvers.f90 decimal.o tipos.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./RiemannSolvers.f90
correctors.o: ./correctors.f90 decimal.o tipos.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./correctors.f90
decimal.o: ./decimal.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./decimal.f90
fluxes.o: ./fluxes.f90 decimal.o tipos.o RiemannSolvers.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluxes.f90
main.o: ./main.f90 decimal.o fluxes.o tipos.o correctors.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./main.f90
plot.o: ./plot.f90 decimal.o tipos.o util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./plot.f90
tipos.o: ./tipos.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./tipos.f90
util.o: ./util.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./util.f90
SRC = ./tipos.f90 ./RiemannSolvers.f90 ./plot.f90 ./decimal.f90 ./correctors.f90 ./fluxes.f90 ./main.f90 ./util.f90
OBJ = tipos.o RiemannSolvers.o plot.o decimal.o correctors.o fluxes.o main.o util.o
clean: neat
	-rm -f .../bin/SWSolver.cppdefs $(OBJ) ../bin/SWSolver
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
../bin/SWSolver: $(OBJ) 
	$(LD) $(OBJ) -o ../bin/SWSolver  $(LDFLAGS)
