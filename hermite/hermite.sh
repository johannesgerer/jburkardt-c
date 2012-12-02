#!/bin/bash
#
cp hermite.h /$HOME/include
#
gcc -c -g -I /$HOME/include hermite.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite.c."
  exit
fi
rm compiler.txt
#
mv hermite.o ~/libc/$ARCH/hermite.o
#
echo "Library installed as ~/libc/$ARCH/hermite.o"
