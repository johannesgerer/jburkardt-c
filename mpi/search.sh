#!/bin/bash
#
mpicc -c search_mpi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling search_mpi.c"
  exit
fi
rm compiler.txt
#
mpicc search_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading search_mpi.o"
  exit
fi
rm search_mpi.o
#
mv a.out search
mpirun -v -np 4 ./search > search_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running search"
  exit
fi
rm search
#
echo "Program output written to search_output.txt"
