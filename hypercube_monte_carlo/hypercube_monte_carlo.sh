#!/bin/bash
#
cp hypercube_monte_carlo.h /$HOME/include
#
gcc -c -g -I /$HOME/include hypercube_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hypercube_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv hypercube_monte_carlo.o ~/libc/$ARCH/hypercube_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/hypercube_monte_carlo.o"
