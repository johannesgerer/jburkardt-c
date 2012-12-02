#!/bin/bash
#
cp unicycle.h /$HOME/include
#
gcc -c -g unicycle.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling unicycle.c"
  exit
fi
rm compiler.txt
#
mv unicycle.o ~/libc/$ARCH/unicycle.o
#
echo "Library installed as ~/libc/$ARCH/unicycle.o"
