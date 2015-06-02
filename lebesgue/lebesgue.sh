#!/bin/bash
#
cp lebesgue.h /$HOME/include
#
gcc -c -g -I/$HOME/include lebesgue.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lebesgue.c"
  exit
fi
rm compiler.txt
#
mv lebesgue.o ~/libc/$ARCH/lebesgue.o
#
echo "Library installed as ~/libc/$ARCH/lebesgue.o"
