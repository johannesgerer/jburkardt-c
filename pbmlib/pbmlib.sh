#!/bin/bash
#
cp pbmlib.h /$HOME/include
#
gcc -c -g pbmlib.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pbmlib.c."
  exit
fi
rm compiler.txt
#
mv pbmlib.o ~/libc/$ARCH/pbmlib.o
#
echo "Library installed as ~/libc/$ARCH/pbmlib.o"
