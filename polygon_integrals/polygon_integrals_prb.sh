#!/bin/bash
#
gcc -c -g -I/$HOME/include polygon_integrals_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_integrals_prb.c"
  exit
fi
rm compiler.txt
#
gcc polygon_integrals_prb.o /$HOME/libc/$ARCH/polygon_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading polygon_integrals_prb.o."
  exit
fi
#
rm polygon_integrals_prb.o
#
mv a.out polygon_integrals_prb
./polygon_integrals_prb > polygon_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running polygon_integrals_prb."
  exit
fi
rm polygon_integrals_prb
#
echo "Program output written to polygon_integrals_prb_output.txt"
