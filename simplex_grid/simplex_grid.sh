#!/bin/bash
#
cp simplex_grid.h /$HOME/include
#
gcc -c -I/$HOME/include simplex_grid.c
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_grid.c"
  exit
fi
#
mv simplex_grid.o ~/libc/$ARCH/simplex_grid.o
#
echo "Library installed as ~/libc/$ARCH/simplex_grid.o"
