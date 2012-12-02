#!/bin/bash
#
cp lagrange_approx_1d.h /$HOME/include
#
gcc -c -g -I/$HOME/include lagrange_approx_1d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lagrange_approx_1d.c"
  exit
fi
rm compiler.txt
#
mv lagrange_approx_1d.o ~/libc/$ARCH/lagrange_approx_1d.o
#
echo "Library installed as ~/libc/$ARCH/lagrange_approx_1d.o"
