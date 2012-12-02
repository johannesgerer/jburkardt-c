#!/bin/bash
#
cp mm_io.h /$HOME/include
#
gcc -c -g mm_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mm_io.c."
  exit
fi
rm compiler.txt
#
mv mm_io.o ~/libc/$ARCH/mm_io.o
#
echo "Library installed as ~/libc/$ARCH/mm_io.o"
