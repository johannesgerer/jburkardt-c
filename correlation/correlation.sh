#!/bin/bash
#
cp correlation.h /$HOME/include
#
gcc -c -g -I/$HOME/include correlation.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling correlation.c"
  exit
fi
rm compiler.txt
#
mv correlation.o ~/libc/$ARCH/correlation.o
#
echo "Library installed as ~/libc/$ARCH/correlation.o"
