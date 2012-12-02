#!/bin/bash
#
cp mgmres.h /$HOME/include
#
gcc -c -g mgmres.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mgmres.c."
  exit
fi
rm compiler.txt
#
mv mgmres.o ~/libc/$ARCH/mgmres.o
#
echo "Library installed as ~/libc/$ARCH/mgmres.o"
