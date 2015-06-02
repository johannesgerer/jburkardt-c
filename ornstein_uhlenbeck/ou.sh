#!/bin/bash
#
cp ou.h /$HOME/include
#
gcc -c -g -I/$HOME/include ou.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ou.c"
  exit
fi
rm compiler.txt
#
mv ou.o ~/libc/$ARCH/ou.o
#
echo "Library installed as ~/libc/$ARCH/ou.o"
