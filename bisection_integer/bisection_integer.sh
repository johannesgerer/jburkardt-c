#!/bin/bash
#
cp bisection_integer.h /$HOME/include
#
gcc -c -g bisection_integer.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bisection_integer.c"
  exit
fi
rm compiler.txt
#
mv bisection_integer.o ~/libc/$ARCH/bisection_integer.o
#
echo "Library installed as ~/libc/$ARCH/bisection_integer.o"
