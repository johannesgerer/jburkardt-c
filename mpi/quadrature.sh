#!/bin/bash
#
mpicc -c quadrature_mpi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling quadrature_mpi.c"
  exit
fi
rm compiler.txt
#
mpicc quadrature_mpi.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading quadrature_mpi.o"
  exit
fi
rm quadrature_mpi.o
#
mv a.out quadrature
mpirun -v -np 4 ./quadrature > quadrature_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running quadrature"
  exit
fi
rm quadrature
#
echo "Program output written to quadrature_output.txt"
