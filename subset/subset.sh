#!/bin/bash
#
cp subset.h /$HOME/include
#
gcc -c -g subset.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subset.c."
  exit
fi
rm compiler.txt
#
mv subset.o ~/libc/$ARCH/subset.o
#
echo "Library installed as ~/libc/$ARCH/subset.o"
