#!/bin/bash
#
cp hermite_polynomial.h /$HOME/include
#
gcc -c -I /$HOME/include hermite_polynomial.c
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_polynomial.c"
  exit
fi
#
mv hermite_polynomial.o ~/libc/$ARCH/hermite_polynomial.o
#
echo "Library installed as ~/libc/$ARCH/hermite_polynomial.o"
