#!/bin/bash
#
gcc -c -g -I/$HOME/include polygon_moments_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_moments_prb.c"
  exit
fi
rm compiler.txt
#
gcc polygon_moments_prb.o /$HOME/libc/$ARCH/polygon_moments.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading polygon_moments_prb.o."
  exit
fi
#
rm polygon_moments_prb.o
#
mv a.out polygon_moments_prb
./polygon_moments_prb > polygon_moments_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running polygon_moments_prb."
  exit
fi
rm polygon_moments_prb
#
echo "Program output written to polygon_moments_prb_output.txt"
