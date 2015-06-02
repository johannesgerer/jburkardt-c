#!/bin/bash
#
mpicc -c -g wave_mpi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wave_mpi.c."
  exit
fi
rm compiler.txt
#
mpicc wave_mpi.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wave_mpi.o."
  exit
fi
#
rm wave_mpi.o
#
mv a.out wave_mpi
mpirun -np 4 ./wave_mpi > wave_local_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wave_mpi"
  exit
fi
rm wave_mpi
#
echo "Program output written to wave_local_output.txt"
