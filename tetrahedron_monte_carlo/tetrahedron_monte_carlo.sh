#!/bin/bash
#
cp tetrahedron_monte_carlo.h /$HOME/include
#
gcc -c -g -I /$HOME/include tetrahedron_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv tetrahedron_monte_carlo.o ~/libc/$ARCH/tetrahedron_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/tetrahedron_monte_carlo.o"
