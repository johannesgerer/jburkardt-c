#!/bin/bash
#
gcc -c -g matmul.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling matmul.c."
  exit
fi
rm compiler.txt
#
gcc matmul.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading matmul.o."
  exit
fi
#
rm matmul.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/matmul
#
echo "Executable installed as ~/binc/$ARCH/matmul"
