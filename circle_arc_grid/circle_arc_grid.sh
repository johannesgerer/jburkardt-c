#!/bin/bash
#
cp circle_arc_grid.h /$HOME/include
#
gcc -c -g circle_arc_grid.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_arc_grid.c."
  exit
fi
rm compiler.txt
#
mv circle_arc_grid.o ~/libc/$ARCH/circle_arc_grid.o
#
echo "Library installed as ~/libc/$ARCH/circle_arc_grid.o"
