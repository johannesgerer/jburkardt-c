#!/bin/bash
#
cp pbmb_io.h /$HOME/include
#
gcc -c -g pbmb_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pbmb_io.c"
  exit
fi
rm compiler.txt
#
mv pbmb_io.o ~/libc/$ARCH/pbmb_io.o
#
echo "Library installed as ~/libc/$ARCH/pbmb_io.o"
