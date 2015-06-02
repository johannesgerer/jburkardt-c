#!/bin/bash
#
cp sphere_monte_carlo.h /$HOME/include
#
gcc -c -g -I /$HOME/include sphere_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv sphere_monte_carlo.o ~/libc/$ARCH/sphere_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/sphere_monte_carlo.o"
