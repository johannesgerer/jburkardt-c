#!/bin/bash
#
cp divdif.h /$HOME/include
#
gcc -c -g -I /$HOME/include divdif.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling divdif.c."
  exit
fi
rm compiler.txt
#
mv divdif.o ~/libc/$ARCH/divdif.o
#
echo "Library installed as ~/libc/$ARCH/divdif.o"
