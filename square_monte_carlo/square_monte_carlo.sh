#!/bin/bash
#
cp square_monte_carlo.h /$HOME/include
#
gcc -c -g -I /$HOME/include square_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling square_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv square_monte_carlo.o ~/libc/$ARCH/square_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/square_monte_carlo.o"
