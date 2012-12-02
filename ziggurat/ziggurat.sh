#!/bin/bash
#
cp ziggurat.h /$HOME/include
#
gcc -c -g ziggurat.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ziggurat.c."
  exit
fi
rm compiler.txt
#
mv ziggurat.o ~/libc/$ARCH/ziggurat.o
#
echo "Library installed as ~/libc/$ARCH/ziggurat.o"
