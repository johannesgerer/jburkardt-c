#!/bin/bash
#
cp i8lib.h /$HOME/include
#
gcc -c -g i8lib.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling i8lib.c."
  exit
fi
rm compiler.txt
#
mv i8lib.o ~/libc/$ARCH/i8lib.o
#
echo "Library installed as ~/libc/$ARCH/i8lib.o"
