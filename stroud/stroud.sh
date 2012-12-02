#!/bin/bash
#
cp stroud.h /$HOME/include
#
gcc -c -g stroud.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stroud.c."
  exit
fi
rm compiler.txt
#
mv stroud.o ~/libc/$ARCH/stroud.o
#
echo "Library installed as ~/libc/$ARCH/stroud.o"
