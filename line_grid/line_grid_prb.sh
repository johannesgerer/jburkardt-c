#!/bin/bash
#
gcc -c -I/$HOME/include line_grid_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling line_grid_prb.c"
  exit
fi
#
gcc -o line_grid_prb line_grid_prb.o /$HOME/libc/$ARCH/line_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading line_grid_prb.o."
  exit
fi
#
rm line_grid_prb.o
#
./line_grid_prb > line_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running line_grid_prb."
  exit
fi
rm line_grid_prb
#
echo "Program output written to line_grid_prb_output.txt"
