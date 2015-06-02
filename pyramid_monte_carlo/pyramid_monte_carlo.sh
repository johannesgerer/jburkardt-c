#!/bin/bash
#
cp pyramid_monte_carlo.h /$HOME/include
#
gcc -c -I/$HOME/include pyramid_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv pyramid_monte_carlo.o ~/libc/$ARCH/pyramid_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/pyramid_monte_carlo.o"
