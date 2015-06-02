#!/bin/bash
#
cp polygon_monte_carlo.h /$HOME/include
#
gcc -c -I/$HOME/include polygon_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv polygon_monte_carlo.o ~/libc/$ARCH/polygon_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/polygon_monte_carlo.o"
