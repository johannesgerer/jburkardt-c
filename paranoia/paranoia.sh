#!/bin/bash
#
gcc -c paranoia.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling paranoia.c"
  exit
fi
rm compiler.txt
#
gcc paranoia.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking paranoia.o"
  exit
fi
#
rm paranoia.o
mv a.out ~/binc/$ARCH/paranoia
#
echo "Executable installed as ~/binc/$ARCH/paranoia"
