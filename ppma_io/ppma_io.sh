#!/bin/bash
#
cp ppma_io.h /$HOME/include
#
gcc -c -g -I /$HOME/include ppma_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ppma_io.c."
  exit
fi
rm compiler.txt
#
mv ppma_io.o ~/libc/$ARCH/ppma_io.o
#
echo "Library installed as ~/libc/$ARCH/ppma_io.o"
