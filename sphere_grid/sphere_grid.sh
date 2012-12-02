#!/bin/bash
#
cp sphere_grid.h /$HOME/include
#
gcc -c -g -I /$HOME/include sphere_grid.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_grid.c."
  exit
fi
rm compiler.txt
#
mv sphere_grid.o ~/libc/$ARCH/sphere_grid.o
#
echo "Library installed as ~/libc/$ARCH/sphere_grid.o"
