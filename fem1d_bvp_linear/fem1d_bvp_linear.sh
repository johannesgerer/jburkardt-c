#!/bin/bash
#
cp fem1d_bvp_linear.h /$HOME/include
#
gcc -c -I /$HOME/include fem1d_bvp_linear.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_bvp_linear.c"
  exit
fi
#
mv fem1d_bvp_linear.o ~/libc/$ARCH/fem1d_bvp_linear.o
#
echo "Library installed as ~/libc/$ARCH/fem1d_bvp_linear.o"
