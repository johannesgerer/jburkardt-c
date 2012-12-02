#!/bin/bash
#
mpicc -c matvec_mpi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling matvec_mpi.c"
  exit
fi
rm compiler.txt
#
mpicc matvec_mpi.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading matvec_mpi.o"
  exit
fi
rm matvec_mpi.o
#
mv a.out matvec
mpirun -v -np 4 ./matvec > matvec_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running matvec"
  exit
fi
rm matvec
#
echo "Program output written to matvec_output.txt"
