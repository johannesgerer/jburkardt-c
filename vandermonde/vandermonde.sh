#!/bin/bash
#
cp vandermonde.h /$HOME/include
#
gcc -c -g -I/$HOME/include vandermonde.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vandermonde.c"
  exit
fi
rm compiler.txt
#
mv vandermonde.o ~/libc/$ARCH/vandermonde.o
#
echo "Library installed as ~/libc/$ARCH/vandermonde.o"
