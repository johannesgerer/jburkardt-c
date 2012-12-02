#!/bin/bash
#
cp c4lib.h /$HOME/include
#
gcc -c -g c4lib.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c4lib.c."
  exit
fi
rm compiler.txt
#
mv c4lib.o ~/libc/$ARCH/c4lib.o
#
echo "Library installed as ~/libc/$ARCH/c4lib.o"
