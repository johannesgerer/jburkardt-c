#!/bin/bash
#
gcc -c -g disk_grid_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling disk_grid_prb.c."
  exit
fi
rm compiler.txt
#
gcc disk_grid_prb.o /$HOME/libc/$ARCH/disk_grid.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading disk_grid_prb.o."
  exit
fi
#
rm disk_grid_prb.o
#
mv a.out disk_grid_prb
./disk_grid_prb > disk_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running disk_grid_prb."
  exit
fi
rm disk_grid_prb
#
echo "Program output written to disk_grid_prb_output.txt"
