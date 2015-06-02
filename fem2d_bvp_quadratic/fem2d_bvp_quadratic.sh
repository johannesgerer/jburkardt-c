#!/bin/bash
#
cp fem2d_bvp_quadratic.h /$HOME/include
#
gcc -c -I /$HOME/include fem2d_bvp_quadratic.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_bvp_quadratic.c"
  exit
fi
#
mv fem2d_bvp_quadratic.o ~/libc/$ARCH/fem2d_bvp_quadratic.o
#
echo "Library installed as ~/libc/$ARCH/fem2d_bvp_quadratic.o"
