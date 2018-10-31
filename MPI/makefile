# fortran compiler
# FC90 = gfortran
# MPI90 = mpif90
FC90 = ifort
MPI90 = mpifort


# optimization level
FLAGS = -O3

# for debugging
# FLAGS = -g -fcheck=all -fbacktrace -Og

# generate warnings and direct to compile.log
WARN = -Wall  > compile.log 2>&1

# add WARN to FLAGS
# FLAGS += $(WARN)

mpiDDC_files = DDC_mpi.f90 mod_DDC_mpi.f90 kdtree2.f90

all: mpiDDC_exe

mpiDDC_exe: $(mpiDDC_files)
	$(MPI90) kdtree2.f90 -c $(FLAGS)
	$(MPI90) mod_DDC_mpi.f90 -c $(FLAGS)
	$(MPI90) $(mpiDDC_files) -o $@ $(FLAGS)

clean:
	rm -r *_exe *.mod *.o *.out *.err dolomite_chem.in myrun* *.dSYM *.inp compile.log locs.txt mass.txt
