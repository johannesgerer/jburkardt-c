#!/bin/bash
#
gcc -c -g fem1d_pmethod.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_pmethod.c"
  exit
fi
rm compiler.txt
#
gcc fem1d_pmethod.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_pmethod.o."
  exit
fi
#
rm fem1d_pmethod.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/fem1d_pmethod
#
echo "Executable installed as ~/binc/$ARCH/fem1d_pmethod"
