#!/bin/bash
#
cp sphere_llq_grid.h /$HOME/include
#
gcc -c -I/$HOME/include sphere_llq_grid.c
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_llq_grid.c"
  exit
fi
#
mv sphere_llq_grid.o ~/libc/$ARCH/sphere_llq_grid.o
#
echo "Library installed as ~/libc/$ARCH/sphere_llq_grid.o"
