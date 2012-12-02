#!/bin/bash
#
cp cycle_floyd.h /$HOME/include
#
gcc -c -g -I /$HOME/include cycle_floyd.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cycle_floyd.c"
  exit
fi
rm compiler.txt
#
mv cycle_floyd.o ~/libc/$ARCH/cycle_floyd.o
#
echo "Library installed as ~/libc/$ARCH/cycle_floyd.o"
