#!/bin/bash
#
cp fem1d_bvp_quadratic.h /$HOME/include
#
gcc -c -I /$HOME/include fem1d_bvp_quadratic.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_bvp_quadratic.c"
  exit
fi
#
mv fem1d_bvp_quadratic.o ~/libc/$ARCH/fem1d_bvp_quadratic.o
#
echo "Library installed as ~/libc/$ARCH/fem1d_bvp_quadratic.o"
