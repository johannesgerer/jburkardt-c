#!/bin/bash
#
cp sftpack.h /$HOME/include
#
gcc -c -g sftpack.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sftpack.c."
  exit
fi
rm compiler.txt
#
mv sftpack.o ~/libc/$ARCH/sftpack.o
#
echo "Library installed as ~/libc/$ARCH/sftpack.o"
