#!/bin/bash
#
cp kronrod.h /$HOME/include
#
gcc -c -g kronrod.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kronrod.c."
  exit
fi
rm compiler.txt
#
mv kronrod.o ~/libc/$ARCH/kronrod.o
#
echo "Library installed as ~/libc/$ARCH/kronrod.o"
