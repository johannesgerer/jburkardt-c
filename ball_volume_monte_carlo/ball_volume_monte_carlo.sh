#!/bin/bash
#
gcc -c -g ball_volume_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ball_volume_monte_carlo.c"
  exit
fi
rm compiler.txt
#
g++ ball_volume_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ball_volume_monte_carlo.o."
  exit
fi
#
rm ball_volume_monte_carlo.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/ball_volume_monte_carlo
#
echo "Executable installed as ~/binc/$ARCH/ball_volume_monte_carlo"
