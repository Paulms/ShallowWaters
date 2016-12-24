# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include template


.DEFAULT:
	-echo $@ does not exist.
all: ../bin/SWSolver
RiemannSolvers.o: ./RiemannSolvers.f90 decimal.o tipos.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./RiemannSolvers.f90
correctors.o: ./correctors.f90 decimal.o tipos.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./correctors.f90
datos.o: ./datos.f90 decimal.o util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./datos.f90
decimal.o: ./decimal.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./decimal.f90
errores.o: ./errores.f90 decimal.o funciones.o quadrature.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./errores.f90
fluxes.o: ./fluxes.f90 decimal.o tipos.o RiemannSolvers.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./fluxes.f90
funciones.o: ./funciones.f90 decimal.o tipos.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./funciones.f90
main.o: ./main.f90 decimal.o fluxes.o tipos.o correctors.o plot.o shallowwatersutils.o datos.o errores.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./main.f90
plot.o: ./plot.f90 decimal.o tipos.o util.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./plot.f90
quadrature.o: ./quadrature.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./quadrature.f90
shallowwatersutils.o: ./shallowwatersutils.f90 decimal.o tipos.o funciones.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./shallowwatersutils.f90
tipos.o: ./tipos.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./tipos.f90
util.o: ./util.f90 decimal.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./util.f90
SRC = ./util.f90 ./RiemannSolvers.f90 ./errores.f90 ./tipos.f90 ./plot.f90 ./fluxes.f90 ./correctors.f90 ./main.f90 ./decimal.f90 ./datos.f90 ./shallowwatersutils.f90 ./quadrature.f90 ./funciones.f90
OBJ = util.o RiemannSolvers.o errores.o tipos.o plot.o fluxes.o correctors.o main.o decimal.o datos.o shallowwatersutils.o quadrature.o funciones.o
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
