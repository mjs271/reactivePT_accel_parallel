#!/bin/bash

make clean 2> /dev/null
make
# # run the program on xx cores using '-np xx' and redirect the error output
# # may need to oversubscribe using the following line
mpirun -np 2 ./mpiDDC_exe 2> a.err
# mpirun --oversubscribe -np 4 ./mpiDDC_exe 2> a.err
# # run the program and redirect all output
# mpirun -np 2 ./mpiDDC_exe > a.out 2> a.err


