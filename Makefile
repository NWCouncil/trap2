FC=mpif90
FFLAGS=-o

trap2: trap2.F90 
	$(FC) $< -o $@ 
	${RM} trap2.o
