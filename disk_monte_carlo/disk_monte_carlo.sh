#!/bin/bash
#
cp disk_monte_carlo.h /$HOME/include
#
gcc -c -g -I/$HOME/include disk_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling disk_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv disk_monte_carlo.o ~/libc/$ARCH/disk_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/disk_monte_carlo.o"
