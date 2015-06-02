#!/bin/bash
#
cp lobatto_polynomial.h /$HOME/include
#
gcc -c -I/$HOME/include lobatto_polynomial.c
if [ $? -ne 0 ]; then
  echo "Errors compiling lobatto_polynomial.c"
  exit
fi
#
mv lobatto_polynomial.o ~/libc/$ARCH/lobatto_polynomial.o
#
echo "Library installed as ~/libc/$ARCH/lobatto_polynomial.o"
