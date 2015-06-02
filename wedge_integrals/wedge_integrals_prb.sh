#!/bin/bash
#
gcc -c -I/$HOME/include wedge_integrals_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_integrals_prb.c"
  exit
fi
#
gcc -o wedge_integrals_prb wedge_integrals_prb.o /$HOME/libc/$ARCH/wedge_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wedge_integrals_prb.o."
  exit
fi
#
rm wedge_integrals_prb.o
#
./wedge_integrals_prb > wedge_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wedge_integrals_prb."
  exit
fi
rm wedge_integrals_prb
#
echo "Program output written to wedge_integrals_prb_output.txt"
