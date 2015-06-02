#!/bin/bash
#
cp triangle_monte_carlo.h /$HOME/include
#
gcc -c -g -I /$HOME/include triangle_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv triangle_monte_carlo.o ~/libc/$ARCH/triangle_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/triangle_monte_carlo.o"
