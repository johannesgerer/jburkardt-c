#!/bin/bash
#
mpicc -c monte_carlo_mpi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling monte_carlo_mpi.c"
  exit
fi
rm compiler.txt
#
mpicc monte_carlo_mpi.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading monte_carlo_mpi.o"
  exit
fi
rm monte_carlo_mpi.o
#
mv a.out monte_carlo
mpirun -v -np 4 ./monte_carlo > monte_carlo_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running monte_carlo"
  exit
fi
rm monte_carlo
#
echo "Program output written to monte_carlo_output.txt"
