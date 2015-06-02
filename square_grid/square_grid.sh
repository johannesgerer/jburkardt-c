#!/bin/bash
#
cp square_grid.h /$HOME/include
#
gcc -c -I/$HOME/include square_grid.c
if [ $? -ne 0 ]; then
  echo "Errors compiling square_grid.c"
  exit
fi
#
mv square_grid.o ~/libc/$ARCH/square_grid.o
#
echo "Library installed as ~/libc/$ARCH/square_grid.o"
