#!/bin/bash
#
cp cube_grid.h /$HOME/include
#
gcc -c -I/$HOME/include cube_grid.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_grid.c"
  exit
fi
#
mv cube_grid.o ~/libc/$ARCH/cube_grid.o
#
echo "Library installed as ~/libc/$ARCH/cube_grid.o"
