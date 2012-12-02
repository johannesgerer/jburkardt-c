#!/bin/bash
#
gcc -c -I/usr/local/dislin life_grid.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling life_grid.c"
  exit
fi
rm compiler.txt
#
gcc life_grid.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading life_grid.o."
  exit
fi
#
rm life_grid.o
#
mv a.out life_grid
./life_grid > life_grid_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running life_grid."
  exit
fi
rm life_grid
#
echo "Program output written to life_grid_output.txt"
