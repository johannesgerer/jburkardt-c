#!/bin/bash
#
cp disk_grid.h /$HOME/include
#
gcc -c -g disk_grid.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling disk_grid.c."
  exit
fi
rm compiler.txt
#
mv disk_grid.o ~/libc/$ARCH/disk_grid.o
#
echo "Library installed as ~/libc/$ARCH/disk_grid.o"
