#!/bin/bash
#
cp linplus.h /$HOME/include
#
gcc -c -g linplus.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linplus.c."
  exit
fi
rm compiler.txt
#
mv linplus.o ~/libc/$ARCH/linplus.o
#
echo "Library installed as ~/libc/$ARCH/linplus.o"
