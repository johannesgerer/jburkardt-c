#!/bin/bash
#
cp wedge_grid.h /$HOME/include
#
gcc -c -I/$HOME/include wedge_grid.c
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_grid.c"
  exit
fi
#
mv wedge_grid.o ~/libc/$ARCH/wedge_grid.o
#
echo "Library installed as ~/libc/$ARCH/wedge_grid.o"
