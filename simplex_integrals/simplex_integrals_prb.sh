#!/bin/bash
#
gcc -c -g -I/$HOME/include simplex_integrals_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_integrals_prb.c"
  exit
fi
rm compiler.txt
#
gcc simplex_integrals_prb.o /$HOME/libc/$ARCH/simplex_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading simplex_integrals_prb.o."
  exit
fi
#
rm simplex_integrals_prb.o
#
mv a.out simplex_integrals_prb
./simplex_integrals_prb > simplex_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running simplex_integrals_prb."
  exit
fi
rm simplex_integrals_prb
#
echo "Program output written to simplex_integrals_prb_output.txt"
