#!/bin/bash
#
gcc -c -g -I/$HOME/include sphere_quad_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_quad_prb.c."
  exit
fi
rm compiler.txt
#
gcc sphere_quad_prb.o /$HOME/libc/$ARCH/sphere_quad.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_quad_prb.o."
  exit
fi
#
rm sphere_quad_prb.o
#
mv a.out sphere_quad_prb
./sphere_quad_prb > sphere_quad_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sphere_quad_prb."
  exit
fi
rm sphere_quad_prb
#
echo "Program output written to sphere_quad_prb_output.txt"
