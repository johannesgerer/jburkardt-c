#!/bin/bash
#
cp components.h /$HOME/include
#
gcc -c -g -I /$HOME/include components.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling components.c"
  exit
fi
rm compiler.txt
#
mv components.o ~/libc/$ARCH/components.o
#
echo "Library installed as ~/libc/$ARCH/components.o"
