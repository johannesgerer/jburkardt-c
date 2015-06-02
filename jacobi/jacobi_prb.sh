#!/bin/bash
#
gcc -c -g -I/$HOME/include jacobi_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling jacobi_prb.c"
  exit
fi
rm compiler.txt
#
gcc jacobi_prb.o /$HOME/libc/$ARCH/jacobi.o ~/libc/$ARCH/gnuplot_i.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading jacobi_prb.o."
  exit
fi
#
rm jacobi_prb.o
#
mv a.out jacobi_prb
./jacobi_prb > jacobi_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running jacobi_prb."
  exit
fi
rm jacobi_prb
#
echo "Program output written to jacobi_prb_output.txt"
