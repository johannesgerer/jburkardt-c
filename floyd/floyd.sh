#!/bin/bash
#
cp floyd.h /$HOME/include
#
gcc -c -g -I /$HOME/include floyd.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling floyd.c."
  exit
fi
rm compiler.txt
#
mv floyd.o ~/libc/$ARCH/floyd.o
#
echo "Library installed as ~/libc/$ARCH/floyd.o"
