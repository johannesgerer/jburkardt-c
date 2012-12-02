#!/bin/bash
#
gcc -c -g ellipsoid_grid_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipsoid_grid_prb.c."
  exit
fi
rm compiler.txt
#
gcc ellipsoid_grid_prb.o /$HOME/libc/$ARCH/ellipsoid_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ellipsoid_grid_prb.o."
  exit
fi
#
rm ellipsoid_grid_prb.o
#
mv a.out ellipsoid_grid_prb
./ellipsoid_grid_prb > ellipsoid_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ellipsoid_grid_prb."
  exit
fi
rm ellipsoid_grid_prb
#
echo "Program output written to ellipsoid_grid_prb_output.txt"
