#!/bin/bash
#
#  Compile
#
gcc -c -g -fopenmp -I/$HOME/include pclinsol.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pclinsol.c"
  exit
fi
rm compiler.txt
#
#  Link and load
#
gcc -fopenmp pclinsol.o -L/$HOME/lib/$ARCH -L/$HOME/libc/$ARCH -lsuperlu_openmp -lm -lblas
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pclinsol.o"
  exit
fi
rm pclinsol.o
mv a.out pclinsol
#
#  Run with 1 processor.
#
export OMP_NUM_THREADS=1
./pclinsol < cg20_cua.txt > pclinsol_1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pclinsol"
  exit
fi
#
#  Run with 4 processors.
#
export OMP_NUM_THREADS=4
./pclinsol < cg20_cua.txt > pclinsol_4_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pclinsol"
  exit
fi
#
rm pclinsol
#
#  Terminate.
#
echo "Program output written to pclinsol_output.txt"
