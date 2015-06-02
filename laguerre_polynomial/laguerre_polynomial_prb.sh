#!/bin/bash
#
gcc -c -g -I/$HOME/include laguerre_polynomial_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laguerre_polynomial_prb.c"
  exit
fi
rm compiler.txt
#
gcc laguerre_polynomial_prb.o /$HOME/libc/$ARCH/laguerre_polynomial.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading laguerre_polynomial_prb.o"
  exit
fi
#
rm laguerre_polynomial_prb.o
#
mv a.out laguerre_polynomial_prb
./laguerre_polynomial_prb > laguerre_polynomial_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running laguerre_polynomial_prb."
  exit
fi
rm laguerre_polynomial_prb
#
echo "Program output written to laguerre_polynomial_prb_output.txt"
