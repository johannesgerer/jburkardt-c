#!/bin/bash
#
cp geometry.h /$HOME/include
#
gcc -c -g geometry.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling geometry.c."
  exit
fi
rm compiler.txt
#
mv geometry.o ~/libc/$ARCH/geometry.o
#
echo "Library installed as ~/libc/$ARCH/geometry.o"
