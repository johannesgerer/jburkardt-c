#!/bin/bash
#
cp cordic.h /$HOME/include
#
gcc -c -g -I /$HOME/include cordic.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cordic.c"
  exit
fi
rm compiler.txt
#
mv cordic.o ~/libc/$ARCH/cordic.o
#
echo "Library installed as ~/libc/$ARCH/cordic.o"
