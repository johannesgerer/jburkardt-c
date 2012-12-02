#!/bin/bash
#
cp eispack.h /$HOME/include
#
gcc -c -g -I/$HOME/include eispack.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling eispack.c"
  exit
fi
rm compiler.txt
#
mv eispack.o ~/libc/$ARCH/eispack.o
#
echo "Library installed as ~/libc/$ARCH/eispack.o"
