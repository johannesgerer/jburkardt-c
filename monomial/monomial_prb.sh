#!/bin/bash
#
gcc -c -I/$HOME/include monomial_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling monomial_prb.c"
  exit
fi
#
gcc monomial_prb.o /$HOME/libc/$ARCH/monomial.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading monomial_prb.o."
  exit
fi
#
rm monomial_prb.o
#
mv a.out monomial_prb
./monomial_prb > monomial_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running monomial_prb."
  exit
fi
rm monomial_prb
#
echo "Program output written to monomial_prb_output.txt"
