#!/bin/bash
#
mpicc -c type_mpi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling type_mpi.c"
  exit
fi
rm compiler.txt
#
mpicc type_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading type_mpi.o"
  exit
fi
rm type_mpi.o
#
mv a.out type
mpirun -v -np 4 ./type > type_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running type"
  exit
fi
rm type
#
echo "Program output written to type_output.txt"
