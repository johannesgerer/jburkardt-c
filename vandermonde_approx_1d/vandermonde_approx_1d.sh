#!/bin/bash
#
cp vandermonde_approx_1d.h /$HOME/include
#
gcc -c -g -I/$HOME/include vandermonde_approx_1d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vandermonde_approx_1d.c"
  exit
fi
rm compiler.txt
#
mv vandermonde_approx_1d.o ~/libc/$ARCH/vandermonde_approx_1d.o
#
echo "Library installed as ~/libc/$ARCH/vandermonde_approx_1d.o"
