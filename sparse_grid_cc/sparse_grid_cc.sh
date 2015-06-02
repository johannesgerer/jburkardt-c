#!/bin/bash
#
cp sparse_grid_cc.h /$HOME/include
#
gcc -c -g -I /$HOME/include sparse_grid_cc.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_cc.c"
  exit
fi
rm compiler.txt
#
mv sparse_grid_cc.o ~/libc/$ARCH/sparse_grid_cc.o
#
echo "Library installed as ~/libc/$ARCH/sparse_grid_cc.o"
