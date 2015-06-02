#!/bin/bash
#
cp blend.h /$HOME/include
#
gcc -c -g -I /$HOME/include blend.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blend.c"
  exit
fi
rm compiler.txt
#
mv blend.o ~/libc/$ARCH/blend.o
#
echo "Library installed as ~/libc/$ARCH/blend.o"
