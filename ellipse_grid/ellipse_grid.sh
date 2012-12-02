#!/bin/bash
#
cp ellipse_grid.h /$HOME/include
#
gcc -c -g ellipse_grid.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipse_grid.c."
  exit
fi
rm compiler.txt
#
mv ellipse_grid.o ~/libc/$ARCH/ellipse_grid.o
#
echo "Library installed as ~/libc/$ARCH/ellipse_grid.o"
