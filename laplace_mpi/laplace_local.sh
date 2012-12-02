#!/bin/bash
#
mpicc -c -g laplace_mpi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laplace_mpi.c."
  exit
fi
rm compiler.txt
#
mpicc laplace_mpi.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading laplace_mpi.o."
  exit
fi
#
rm laplace_mpi.o
#
mv a.out laplace
mpirun -np 4 ./laplace > laplace_local_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running laplace"
  exit
fi
rm laplace
#
echo "Program output written to laplace_local_output.txt"
