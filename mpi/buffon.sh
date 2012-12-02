#!/bin/bash
#
mpicc -c buffon_mpi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling buffon_mpi.c"
  exit
fi
rm compiler.txt
#
mpicc buffon_mpi.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading buffon_mpi.o"
  exit
fi
rm buffon_mpi.o
#
mv a.out buffon
mpirun -v -np 4 ./buffon > buffon_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running buffon"
  exit
fi
rm buffon
#
echo "Program output written to buffon_output.txt"
