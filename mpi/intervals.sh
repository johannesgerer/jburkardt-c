#!/bin/bash
#
mpicc -c intervals_mpi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling intervals_mpi.c"
  exit
fi
rm compiler.txt
#
mpicc intervals_mpi.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading intervals_mpi.o"
  exit
fi
rm intervals_mpi.o
#
mv a.out intervals
mpirun -v -np 4 ./intervals > intervals_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running intervals"
  exit
fi
rm intervals
#
echo "Program output written to intervals_output.txt"
