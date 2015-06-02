#!/bin/bash
#
cp jacobi_polynomial.h /$HOME/include
#
gcc -c -g -I /$HOME/include jacobi_polynomial.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling jacobi_polynomial.c"
  exit
fi
rm compiler.txt
#
mv jacobi_polynomial.o ~/libc/$ARCH/jacobi_polynomial.o
#
echo "Library installed as ~/libc/$ARCH/jacobi_polynomial.o"
