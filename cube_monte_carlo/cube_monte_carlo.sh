#!/bin/bash
#
cp cube_monte_carlo.h /$HOME/include
#
gcc -c -g -I /$HOME/include cube_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv cube_monte_carlo.o ~/libc/$ARCH/cube_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/cube_monte_carlo.o"
