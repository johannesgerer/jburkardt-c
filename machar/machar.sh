#!/bin/bash
#
cp machar.h /$HOME/include
#
gcc -c -g machar.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling machar.c."
  exit
fi
rm compiler.txt
#
mv machar.o ~/libc/$ARCH/machar.o
#
echo "Library installed as ~/libc/$ARCH/machar.o"
