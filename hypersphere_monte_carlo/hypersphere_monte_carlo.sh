#!/bin/bash
#
cp hypersphere_monte_carlo.h /$HOME/include
#
gcc -c -g -I /$HOME/include hypersphere_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hypersphere_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv hypersphere_monte_carlo.o ~/libc/$ARCH/hypersphere_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/hypersphere_monte_carlo.o"
