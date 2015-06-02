#!/bin/bash
#
gcc -c legendre_exactness.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling legendre_exactness.c"
  exit
fi
rm compiler.txt
#
gcc legendre_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading legendre_exactness.o"
  exit
fi
rm legendre_exactness.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/legendre_exactness
#
echo "Executable installed as ~/binc/$ARCH/legendre_exactness"
