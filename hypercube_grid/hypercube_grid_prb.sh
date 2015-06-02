#!/bin/bash
#
gcc -c -I/$HOME/include hypercube_grid_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling hypercube_grid_prb.c"
  exit
fi
#
gcc -o hypercube_grid_prb hypercube_grid_prb.o /$HOME/libc/$ARCH/hypercube_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hypercube_grid_prb.o."
  exit
fi
#
rm hypercube_grid_prb.o
#
./hypercube_grid_prb > hypercube_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hypercube_grid_prb."
  exit
fi
rm hypercube_grid_prb
#
echo "Program output written to hypercube_grid_prb_output.txt"
