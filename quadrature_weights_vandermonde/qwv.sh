#!/bin/bash
#
cp qwv.h /$HOME/include
#
gcc -c -g -I/$HOME/include qwv.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qwv.c"
  exit
fi
rm compiler.txt
#
mv qwv.o ~/libc/$ARCH/qwv.o
#
echo "Library installed as ~/libc/$ARCH/qwv.o"
