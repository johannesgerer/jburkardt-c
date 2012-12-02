#!/bin/bash
#
mpicc -c bones_mpi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling bones_mpi.c"
  exit
fi
rm compiler.txt
#
mpicc bones_mpi.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading bones_mpi.o"
  exit
fi
rm bones_mpi.o
#
mv a.out bones
mpirun -v -np 4 ./bones > bones_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running bones"
  exit
fi
rm bones
#
echo "Program output written to bones_output.txt"
