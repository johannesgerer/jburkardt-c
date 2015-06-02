#!/bin/bash
#
cp quadmom.h /$HOME/include
#
gcc -c -g -I/$HOME/include quadmom.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadmom.c"
  exit
fi
rm compiler.txt
#
mv quadmom.o ~/libc/$ARCH/quadmom.o
#
echo "Library installed as ~/libc/$ARCH/quadmom.o"
