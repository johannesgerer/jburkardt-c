#!/bin/bash
#
mpicc quad_mpi.c
#
if [ $? -ne 0 ]; then
  echo "Errors compiling quad_mpi.c"
  exit
fi
#
#  Rename the executable.
#
mv a.out quad
#
#  Ask MPI to use 8 processes to run your program.
#
mpirun -np 8 ./quad > quad_local_output.txt
#
#  Clean up.
#
rm quad
#
echo "Program output written to quad_local_output.txt"

