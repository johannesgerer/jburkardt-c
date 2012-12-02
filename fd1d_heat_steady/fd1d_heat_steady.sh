#!/bin/bash
#
cp fd1d_heat_steady.h /$HOME/include
#
gcc -c -g -I /$HOME/include fd1d_heat_steady.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_heat_steady.c."
  exit
fi
rm compiler.txt
#
mv fd1d_heat_steady.o ~/libc/$ARCH/fd1d_heat_steady.o
#
echo "Library installed as ~/libc/$ARCH/fd1d_heat_steady.o"
