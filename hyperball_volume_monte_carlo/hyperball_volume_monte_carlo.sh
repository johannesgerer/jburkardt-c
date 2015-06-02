#!/bin/bash
#
gcc -c -g hyperball_volume_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hyperball_volume_monte_carlo.c"
  exit
fi
rm compiler.txt
#
g++ hyperball_volume_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hyperball_volume_monte_carlo.o."
  exit
fi
#
rm hyperball_volume_monte_carlo.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/hyperball_volume_monte_carlo
#
echo "Executable installed as ~/binc/$ARCH/hyperball_volume_monte_carlo"
