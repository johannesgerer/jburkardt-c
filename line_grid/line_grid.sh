#!/bin/bash
#
cp line_grid.h /$HOME/include
#
gcc -c -I/$HOME/include line_grid.c
if [ $? -ne 0 ]; then
  echo "Errors compiling line_grid.c"
  exit
fi
#
mv line_grid.o ~/libc/$ARCH/line_grid.o
#
echo "Library installed as ~/libc/$ARCH/line_grid.o"
