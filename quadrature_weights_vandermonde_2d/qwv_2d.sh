#!/bin/bash
#
cp qwv_2d.h /$HOME/include
#
gcc -c -g -I/$HOME/include qwv_2d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qwv_2d.c"
  exit
fi
rm compiler.txt
#
mv qwv_2d.o ~/libc/$ARCH/qwv_2d.o
#
echo "Library installed as ~/libc/$ARCH/qwv_2d.o"
