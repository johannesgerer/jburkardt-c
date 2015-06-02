#!/bin/bash
#
cp polygon_grid.h /$HOME/include
#
gcc -c -I/$HOME/include polygon_grid.c
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_grid.c"
  exit
fi
#
mv polygon_grid.o ~/libc/$ARCH/polygon_grid.o
#
echo "Library installed as ~/libc/$ARCH/polygon_grid.o"
