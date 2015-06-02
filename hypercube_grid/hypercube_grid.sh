#!/bin/bash
#
cp hypercube_grid.h /$HOME/include
#
gcc -c -I/$HOME/include hypercube_grid.c
if [ $? -ne 0 ]; then
  echo "Errors compiling hypercube_grid.c"
  exit
fi
#
mv hypercube_grid.o ~/libc/$ARCH/hypercube_grid.o
#
echo "Library installed as ~/libc/$ARCH/hypercube_grid.o"
