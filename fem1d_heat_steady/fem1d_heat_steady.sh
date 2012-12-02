#!/bin/bash
#
cp fem1d_heat_steady.h /$HOME/include
#
gcc -c -g -I /$HOME/include fem1d_heat_steady.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_heat_steady.c."
  exit
fi
rm compiler.txt
#
mv fem1d_heat_steady.o ~/libc/$ARCH/fem1d_heat_steady.o
#
echo "Library installed as ~/libc/$ARCH/fem1d_heat_steady.o"
