#!/bin/bash
#
gcc -c -g nsasm.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nsasm.c."
  exit
fi
rm compiler.txt
#
mv nsasm.o ~/libc/$ARCH/nsasm.o
#
echo "Library installed as ~/libc/$ARCH/nsasm.o"
