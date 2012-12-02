#!/bin/bash
#
cp latinize.h /$HOME/include
#
gcc -c -g -I /$HOME/include latinize.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latinize.c"
  exit
fi
rm compiler.txt
#
mv latinize.o ~/libc/$ARCH/latinize.o
#
echo "Library installed as ~/libc/$ARCH/latinize.o"
