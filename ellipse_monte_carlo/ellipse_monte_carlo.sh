#!/bin/bash
#
cp ellipse_monte_carlo.h /$HOME/include
#
gcc -c -I/$HOME/include ellipse_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipse_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv ellipse_monte_carlo.o ~/libc/$ARCH/ellipse_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/ellipse_monte_carlo.o"
