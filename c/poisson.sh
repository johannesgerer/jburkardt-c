#!/bin/bash
#
gcc -c poisson.c
if [ $? -ne 0 ]; then
  echo "Errors while compiling poisson.c"
  exit
fi
#
gcc poisson.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading poisson.o"
  exit
fi
rm poisson.o
#
mv a.out poisson
./poisson > poisson_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running poisson"
  exit
fi
rm poisson
#
echo "Program output written to poisson_output.txt"
