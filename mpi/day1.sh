#!/bin/bash
#
mpicc -c day1_mpi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling day1_mpi.c"
  exit
fi
rm compiler.txt
#
mpicc day1_mpi.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading day1_mpi.o"
  exit
fi
rm day1_mpi.o
#
mv a.out day1
mpirun -v -np 4 ./day1 > day1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running day1"
  exit
fi
rm day1
#
echo "Program output written to day1_output.txt"
