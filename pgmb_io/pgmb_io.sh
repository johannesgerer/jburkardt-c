#!/bin/bash
#
cp pgmb_io.h /$HOME/include
#
gcc -c -g -I /$HOME/include pgmb_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pgmb_io.c."
  exit
fi
rm compiler.txt
#
mv pgmb_io.o ~/libc/$ARCH/pgmb_io.o
#
echo "Library installed as ~/libc/$ARCH/pgmb_io.o"
