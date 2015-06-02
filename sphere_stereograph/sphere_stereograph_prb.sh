#!/bin/bash
#
gcc -c -g -I/$HOME/include sphere_stereograph_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_stereograph_prb.c"
  exit
fi
rm compiler.txt
#
gcc sphere_stereograph_prb.o /$HOME/libc/$ARCH/sphere_stereograph.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_stereograph_prb.o."
  exit
fi
#
rm sphere_stereograph_prb.o
#
mv a.out sphere_stereograph_prb
./sphere_stereograph_prb > sphere_stereograph_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sphere_stereograph_prb."
  exit
fi
rm sphere_stereograph_prb
#
echo "Program output written to sphere_stereograph_prb_output.txt"
