#!/bin/bash
#
mpicc satisfy_mpi.c
#
if [ $? -ne 0 ]; then
  echo "Errors compiling satisfy_mpi.c"
  exit
fi
#
#  Rename the executable.
#
mv a.out satisfy
#
#  Ask MPI to use 4 processes to run your program.
#
mpirun -np 4 ./satisfy > satisfy_local_output.txt
#
#  Clean up.
#
rm satisfy
#
echo "Program output written to satisfy_local_output.txt"

