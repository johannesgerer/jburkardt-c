#!/bin/bash
#
cp vandermonde_interp_2d.h /$HOME/include
#
gcc -c -g -I /$HOME/include vandermonde_interp_2d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vandermonde_interp_2d.c"
  exit
fi
rm compiler.txt
#
mv vandermonde_interp_2d.o ~/libc/$ARCH/vandermonde_interp_2d.o
#
echo "Library installed as ~/libc/$ARCH/vandermonde_interp_2d.o"
