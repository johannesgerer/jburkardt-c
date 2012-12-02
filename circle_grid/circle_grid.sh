#!/bin/bash
#
cp circle_grid.h /$HOME/include
#
gcc -c -g circle_grid.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_grid.c."
  exit
fi
rm compiler.txt
#
mv circle_grid.o ~/libc/$ARCH/circle_grid.o
#
echo "Library installed as ~/libc/$ARCH/circle_grid.o"
