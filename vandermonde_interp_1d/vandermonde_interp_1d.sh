#!/bin/bash
#
cp vandermonde_interp_1d.h /$HOME/include
#
gcc -c -g -I/$HOME/include vandermonde_interp_1d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vandermonde_interp_1d.c"
  exit
fi
rm compiler.txt
#
mv vandermonde_interp_1d.o ~/libc/$ARCH/vandermonde_interp_1d.o
#
echo "Library installed as ~/libc/$ARCH/vandermonde_interp_1d.o"
