#!/bin/bash
#
cp brent.h /$HOME/include
#
gcc -c -g brent.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling brent.c."
  exit
fi
rm compiler.txt
#
mv brent.o ~/libc/$ARCH/brent.o
#
echo "Library installed as ~/libc/$ARCH/brent.o"
