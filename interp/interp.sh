#!/bin/bash
#
cp interp.h /$HOME/include
#
gcc -c -g -I/$HOME/include interp.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling interp.c"
  exit
fi
rm compiler.txt
#
mv interp.o ~/libc/$ARCH/interp.o
#
echo "Library installed as ~/libc/$ARCH/interp.o"
