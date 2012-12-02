#!/bin/bash
#
cp bio.h /$HOME/include
#
gcc -c -g bio.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bio.c."
  exit
fi
rm compiler.txt
#
mv bio.o ~/libc/$ARCH/bio.o
#
echo "Library installed as ~/libc/$ARCH/bio.o"
