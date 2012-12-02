#!/bin/bash
#
cp sde.h /$HOME/include
#
gcc -c -g -I/$HOME/include sde.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sde.c"
  exit
fi
rm compiler.txt
#
mv sde.o ~/libc/$ARCH/sde.o
#
echo "Library installed as ~/libc/$ARCH/sde.o"
