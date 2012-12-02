#!/bin/bash
#
cp ellipsoid_grid.h /$HOME/include
#
gcc -c -g ellipsoid_grid.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipsoid_grid.c."
  exit
fi
rm compiler.txt
#
mv ellipsoid_grid.o ~/libc/$ARCH/ellipsoid_grid.o
#
echo "Library installed as ~/libc/$ARCH/ellipsoid_grid.o"
