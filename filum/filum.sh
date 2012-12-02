#!/bin/bash
#
cp filum.h /$HOME/include
#
gcc -c -g filum.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling filum.c."
  exit
fi
rm compiler.txt
#
mv filum.o ~/libc/$ARCH/filum.o
#
echo "Library installed as ~/libc/$ARCH/filum.o"
