#!/bin/bash
#
cp sphere_fibonacci_grid.h /$HOME/include
#
gcc -c -I/$HOME/include sphere_fibonacci_grid.c
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_fibonacci_grid.c"
  exit
fi
#
mv sphere_fibonacci_grid.o ~/libc/$ARCH/sphere_fibonacci_grid.o
#
echo "Library installed as ~/libc/$ARCH/sphere_fibonacci_grid.o"
