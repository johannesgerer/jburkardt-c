#!/bin/bash
#
cp rk4.h /$HOME/include
#
gcc -c -g -I /$HOME/include rk4.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rk4.c"
  exit
fi
rm compiler.txt
#
mv rk4.o ~/libc/$ARCH/rk4.o
#
echo "Library installed as ~/libc/$ARCH/rk4.o"
