#!/bin/bash
#
gcc -c -g triangulation_triangle_neighbors.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_triangle_neighbors.c."
  exit
fi
rm compiler.txt
#
gcc triangulation_triangle_neighbors.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_triangle_neighbors.o."
  exit
fi
#
rm triangulation_triangle_neighbors.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/triangulation_triangle_neighbors
#
echo "Executable installed as ~/binc/$ARCH/triangulation_triangle_neighbors"
