#!/bin/bash
#
cp simplex_monte_carlo.h /$HOME/include
#
gcc -c -g -I /$HOME/include simplex_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv simplex_monte_carlo.o ~/libc/$ARCH/simplex_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/simplex_monte_carlo.o"
