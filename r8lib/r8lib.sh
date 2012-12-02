#!/bin/bash
#
cp r8lib.h /$HOME/include
#
gcc -c -g -I/$HOME/include r8lib.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling r8lib.c"
  exit
fi
rm compiler.txt
#
mv r8lib.o ~/libc/$ARCH/r8lib.o
#
echo "Library installed as ~/libc/$ARCH/r8lib.o"
