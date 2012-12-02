#!/bin/bash
#
cp timestamp.h /$HOME/include
#
gcc -c -g timestamp.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling timestamp.c."
  exit
fi
rm compiler.txt
#
mv timestamp.o ~/libc/$ARCH/timestamp.o
#
echo "Library installed as ~/libc/$ARCH/timestamp.o"
