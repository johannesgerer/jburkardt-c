#!/bin/bash
#
gcc -c -I/$HOME/include square_grid_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling square_grid_prb.c"
  exit
fi
#
gcc -o square_grid_prb square_grid_prb.o /$HOME/libc/$ARCH/square_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading square_grid_prb.o."
  exit
fi
#
rm square_grid_prb.o
#
./square_grid_prb > square_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running square_grid_prb."
  exit
fi
rm square_grid_prb
#
echo "Program output written to square_grid_prb_output.txt"
