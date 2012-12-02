#!/bin/bash
#
cp c8lib.h /$HOME/include
#
gcc -c -g c8lib.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c8lib.c."
  exit
fi
rm compiler.txt
#
mv c8lib.o ~/libc/$ARCH/c8lib.o
#
echo "Library installed as ~/libc/$ARCH/c8lib.o"
