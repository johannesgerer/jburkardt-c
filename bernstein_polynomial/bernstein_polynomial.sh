#!/bin/bash
#
cp bernstein_polynomial.h /$HOME/include
#
gcc -c -g -I /$HOME/include bernstein_polynomial.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bernstein_polynomial.c"
  exit
fi
rm compiler.txt
#
mv bernstein_polynomial.o ~/libc/$ARCH/bernstein_polynomial.o
#
echo "Library installed as ~/libc/$ARCH/bernstein_polynomial.o"
