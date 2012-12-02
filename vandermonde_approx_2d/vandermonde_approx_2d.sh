#!/bin/bash
#
cp vandermonde_approx_2d.h /$HOME/include
#
gcc -c -g -I/$HOME/include vandermonde_approx_2d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vandermonde_approx_2d.c"
  exit
fi
rm compiler.txt
#
mv vandermonde_approx_2d.o ~/libc/$ARCH/vandermonde_approx_2d.o
#
echo "Library installed as ~/libc/$ARCH/vandermonde_approx_2d.o"
