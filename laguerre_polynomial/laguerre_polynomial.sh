#!/bin/bash
#
cp laguerre_polynomial.h /$HOME/include
#
gcc -c -g -I /$HOME/include laguerre_polynomial.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laguerre_polynomial.c"
  exit
fi
rm compiler.txt
#
mv laguerre_polynomial.o ~/libc/$ARCH/laguerre_polynomial.o
#
echo "Library installed as ~/libc/$ARCH/laguerre_polynomial.o"
