#!/bin/bash
#
gcc -c -g fem1d_nonlinear.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_nonlinear.c"
  exit
fi
rm compiler.txt
#
gcc fem1d_nonlinear.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_nonlinear.o."
  exit
fi
#
rm fem1d_nonlinear.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/fem1d_nonlinear
#
echo "Executable installed as ~/binc/$ARCH/fem1d_nonlinear"
