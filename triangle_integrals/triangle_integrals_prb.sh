#!/bin/bash
#
gcc -c -I/$HOME/include triangle_integrals_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_integrals_prb.c"
  exit
fi
#
gcc -o triangle_integrals_prb triangle_integrals_prb.o /$HOME/libc/$ARCH/triangle_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_integrals_prb.o."
  exit
fi
#
rm triangle_integrals_prb.o
#
./triangle_integrals_prb > triangle_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running triangle_integrals_prb."
  exit
fi
rm triangle_integrals_prb
#
echo "Program output written to triangle_integrals_prb_output.txt"
