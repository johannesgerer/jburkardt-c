#!/bin/bash
#
cp tetrahedron_grid.h /$HOME/include
#
gcc -c -g tetrahedron_grid.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_grid.c."
  exit
fi
rm compiler.txt
#
mv tetrahedron_grid.o ~/libc/$ARCH/tetrahedron_grid.o
#
echo "Library installed as ~/libc/$ARCH/tetrahedron_grid.o"
