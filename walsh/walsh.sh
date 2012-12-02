#!/bin/bash
#
cp walsh.h /$HOME/include
#
gcc -c -g -I /$HOME/include walsh.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling walsh.c."
  exit
fi
rm compiler.txt
#
mv walsh.o ~/libc/$ARCH/walsh.o
#
echo "Library installed as ~/libc/$ARCH/walsh.o"
