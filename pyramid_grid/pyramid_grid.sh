#!/bin/bash
#
cp pyramid_grid.h /$HOME/include
#
gcc -c -I/$HOME/include pyramid_grid.c
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_grid.c"
  exit
fi
#
mv pyramid_grid.o ~/libc/$ARCH/pyramid_grid.o
#
echo "Library installed as ~/libc/$ARCH/pyramid_grid.o"
