#!/bin/bash
#
gcc -c -g fem1d_adaptive.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_adaptive.c"
  exit
fi
rm compiler.txt
#
gcc fem1d_adaptive.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_adaptive.o."
  exit
fi
#
rm fem1d_adaptive.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/fem1d_adaptive
#
echo "Executable installed as ~/binc/$ARCH/fem1d_adaptive"
