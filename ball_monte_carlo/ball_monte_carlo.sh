#!/bin/bash
#
cp ball_monte_carlo.h /$HOME/include
#
gcc -c -g -I/$HOME/include ball_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ball_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv ball_monte_carlo.o ~/libc/$ARCH/ball_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/ball_monte_carlo.o"
