#!/bin/bash
#
cp wedge_monte_carlo.h /$HOME/include
#
gcc -c -I/$HOME/include wedge_monte_carlo.c
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_monte_carlo.c"
  exit
fi
#
mv wedge_monte_carlo.o ~/libc/$ARCH/wedge_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/wedge_monte_carlo.o"
