#!/bin/bash
#
cp circle_monte_carlo.h /$HOME/include
#
gcc -c -g -I /$HOME/include circle_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv circle_monte_carlo.o ~/libc/$ARCH/circle_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/circle_monte_carlo.o"
