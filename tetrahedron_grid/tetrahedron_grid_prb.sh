#!/bin/bash
#
gcc -c -g tetrahedron_grid_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_grid_prb.c."
  exit
fi
rm compiler.txt
#
gcc tetrahedron_grid_prb.o /$HOME/libc/$ARCH/tetrahedron_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tetrahedron_grid_prb.o."
  exit
fi
#
rm tetrahedron_grid_prb.o
#
mv a.out tetrahedron_grid_prb
./tetrahedron_grid_prb > tetrahedron_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tetrahedron_grid_prb."
  exit
fi
rm tetrahedron_grid_prb
#
echo "Program output written to tetrahedron_grid_prb_output.txt"
