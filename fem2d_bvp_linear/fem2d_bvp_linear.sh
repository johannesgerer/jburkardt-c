#!/bin/bash
#
cp fem2d_bvp_linear.h /$HOME/include
#
gcc -c -I /$HOME/include fem2d_bvp_linear.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_bvp_linear.c"
  exit
fi
#
mv fem2d_bvp_linear.o ~/libc/$ARCH/fem2d_bvp_linear.o
#
echo "Library installed as ~/libc/$ARCH/fem2d_bvp_linear.o"
