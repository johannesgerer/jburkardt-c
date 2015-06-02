#!/bin/bash
#
cp line_monte_carlo.h /$HOME/include
#
gcc -c -g -I /$HOME/include line_monte_carlo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling line_monte_carlo.c"
  exit
fi
rm compiler.txt
#
mv line_monte_carlo.o ~/libc/$ARCH/line_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/line_monte_carlo.o"
