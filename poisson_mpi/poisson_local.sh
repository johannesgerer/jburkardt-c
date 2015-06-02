#!/bin/bash
#
mpicc poisson_mpi.c
#
if [ $? -ne 0 ]; then
  echo "Errors compiling poisson_mpi.c"
  exit
fi
#
#  Rename the executable.
#
mv a.out poisson
#
#  Ask MPI to use 8 processes to run the program.
#
mpirun -np 8 ./poisson > poisson_local_output.txt
#
#  Clean up.
#
rm poisson
#
echo "Program output written to poisson_local_output.txt"

