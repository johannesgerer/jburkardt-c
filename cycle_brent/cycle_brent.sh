#!/bin/bash
#
cp cycle_brent.h /$HOME/include
#
gcc -c -g -I /$HOME/include cycle_brent.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cycle_brent.c"
  exit
fi
rm compiler.txt
#
mv cycle_brent.o ~/libc/$ARCH/cycle_brent.o
#
echo "Library installed as ~/libc/$ARCH/cycle_brent.o"
