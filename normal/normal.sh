#!/bin/bash
#
cp normal.h /$HOME/include
#
gcc -c -g normal.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling normal.c."
  exit
fi
rm compiler.txt
#
mv normal.o ~/libc/$ARCH/normal.o
#
echo "Library installed as ~/libc/$ARCH/normal.o"
