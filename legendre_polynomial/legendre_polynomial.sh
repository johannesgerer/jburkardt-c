#!/bin/bash
#
cp legendre_polynomial.h /$HOME/include
#
gcc -c -I /$HOME/include legendre_polynomial.c
if [ $? -ne 0 ]; then
  echo "Errors compiling legendre_polynomial.c"
  exit
fi
#
mv legendre_polynomial.o ~/libc/$ARCH/legendre_polynomial.o
#
echo "Library installed as ~/libc/$ARCH/legendre_polynomial.o"
