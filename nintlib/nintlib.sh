#!/bin/bash
#
cp nintlib.h /$HOME/include
#
gcc -c -g -I /$HOME/include nintlib.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nintlib.c"
  exit
fi
rm compiler.txt
#
mv nintlib.o ~/libc/$ARCH/nintlib.o
#
echo "Library installed as ~/libc/$ARCH/nintlib.o"
