#!/bin/bash
#
cp sphere_llt_grid.h /$HOME/include
#
gcc -c -I/$HOME/include sphere_llt_grid.c
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_llt_grid.c"
  exit
fi
#
mv sphere_llt_grid.o ~/libc/$ARCH/sphere_llt_grid.o
#
echo "Library installed as ~/libc/$ARCH/sphere_llt_grid.o"
