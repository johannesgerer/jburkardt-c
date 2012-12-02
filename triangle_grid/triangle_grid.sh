#!/bin/bash
#
cp triangle_grid.h /$HOME/include
#
gcc -c -g -I /$HOME/include triangle_grid.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_grid.c."
  exit
fi
rm compiler.txt
#
mv triangle_grid.o ~/libc/$ARCH/triangle_grid.o
#
echo "Library installed as ~/libc/$ARCH/triangle_grid.o"
