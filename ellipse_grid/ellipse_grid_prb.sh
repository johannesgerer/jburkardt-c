#!/bin/bash
#
gcc -c -g ellipse_grid_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipse_grid_prb.c."
  exit
fi
rm compiler.txt
#
gcc ellipse_grid_prb.o /$HOME/libc/$ARCH/ellipse_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ellipse_grid_prb.o."
  exit
fi
#
rm ellipse_grid_prb.o
#
mv a.out ellipse_grid_prb
./ellipse_grid_prb > ellipse_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ellipse_grid_prb."
  exit
fi
rm ellipse_grid_prb
#
echo "Program output written to ellipse_grid_prb_output.txt"
