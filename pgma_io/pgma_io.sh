#!/bin/bash
#
cp pgma_io.h /$HOME/include
#
gcc -c -g pgma_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pgma_io.c."
  exit
fi
rm compiler.txt
#
mv pgma_io.o ~/libc/$ARCH/pgma_io.o
#
echo "Library installed as ~/libc/$ARCH/pgma_io.o"
