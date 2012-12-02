#!/bin/bash
#
cp i4lib.h /$HOME/include
#
gcc -c -g i4lib.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling i4lib.c."
  exit
fi
rm compiler.txt
#
mv i4lib.o ~/libc/$ARCH/i4lib.o
#
echo "Library installed as ~/libc/$ARCH/i4lib.o"
