#!/bin/bash
#
cp hyperball_monte_carlo.h /$HOME/include
#
gcc -c -g -I/$HOME/include hyperball_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hyperball_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv hyperball_monte_carlo.o ~/libc/$ARCH/hyperball_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/hyperball_monte_carlo.o"
