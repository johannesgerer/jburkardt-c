#!/bin/bash
#
cp cyclic_reduction.h /$HOME/include
#
gcc -c -g cyclic_reduction.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cyclic_reduction.c."
  exit
fi
rm compiler.txt
#
mv cyclic_reduction.o ~/libc/$ARCH/cyclic_reduction.o
#
echo "Library installed as ~/libc/$ARCH/cyclic_reduction.o"
