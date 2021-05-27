FC = gfortran
FFLAGS = -o
CFLAGS = -c -fopenmp
GFLAGS = -pg -g -fopenmp
#CFLAGS = -c -fbounds-check
#GFLAGS = -pg -g -fbounds-check 
%.o: %.mod

all: calcphis

calcphis: procedures.o param.o hamiltonians.o sc_tensor.o calcphis.f90
	$(FC) $(FFLAGS) calcphis procedures.o param.o hamiltonians.o sc_tensor.o calcphis.f90 $(GFLAGS)

sc_tensor.o: sc_tensor.f90 
	$(FC) $(CFLAGS) sc_tensor.f90

hamiltonians.o: hamiltonians.f90
	$(FC) $(CFLAGS) hamiltonians.f90

param.o: param.f90 
	$(FC) $(CFLAGS) param.f90

procedures.o: procedures.f90
	$(FC) $(CFLAGS) procedures.f90

clean:
	rm *o calcphis

