#!/bin/bash
#
gcc -c -g bisection_integer_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bisection_integer_prb.c"
  exit
fi
rm compiler.txt
#
gcc bisection_integer_prb.o /$HOME/libc/$ARCH/bisection_integer.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bisection_integer_prb.o."
  exit
fi
#
rm bisection_integer_prb.o
#
mv a.out bisection_integer_prb
./bisection_integer_prb > bisection_integer_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bisection_integer_prb."
  exit
fi
rm bisection_integer_prb
#
echo "Program output written to bisection_integer_prb_output.txt"
