#!/bin/bash
#
cp triangle_io.h /$HOME/include
#
gcc -c -g triangle_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_io.c."
  exit
fi
rm compiler.txt
#
mv triangle_io.o ~/libc/$ARCH/triangle_io.o
#
echo "Library installed as ~/libc/$ARCH/triangle_io.o"
