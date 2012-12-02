#!/bin/bash
#
cp pbma_io.h /$HOME/include
#
gcc -c -g pbma_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pbma_io.c."
  exit
fi
rm compiler.txt
#
mv pbma_io.o ~/libc/$ARCH/pbma_io.o
#
echo "Library installed as ~/libc/$ARCH/pbma_io.o"
