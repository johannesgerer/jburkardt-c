#!/bin/bash
#
cp sparse_grid_hw.h /$HOME/include
#
gcc -c -g -I/$HOME/include sparse_grid_hw.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_hw.c"
  exit
fi
rm compiler.txt
#
mv sparse_grid_hw.o ~/libc/$ARCH/sparse_grid_hw.o
#
echo "Library installed as ~/libc/$ARCH/sparse_grid_hw.o"
