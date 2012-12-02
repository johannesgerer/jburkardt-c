#!/bin/bash
#
mpicc -c sum_mpi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling sum_mpi.c"
  exit
fi
rm compiler.txt
#
mpicc sum_mpi.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading sum_mpi.o"
  exit
fi
rm sum_mpi.o
#
mv a.out sum
mpirun -v -np 4 ./sum > sum_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running sum"
  exit
fi
rm sum
#
echo "Program output written to sum_output.txt"
