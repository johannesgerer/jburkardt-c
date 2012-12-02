#!/bin/bash
#
cp box_behnken.h /$HOME/include
#
gcc -c -g -I /$HOME/include box_behnken.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling box_behnken.c"
  exit
fi
rm compiler.txt
#
mv box_behnken.o ~/libc/$ARCH/box_behnken.o
#
echo "Library installed as ~/libc/$ARCH/box_behnken.o"
