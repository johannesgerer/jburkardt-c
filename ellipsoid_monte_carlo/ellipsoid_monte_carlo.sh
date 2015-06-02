#!/bin/bash
#
cp ellipsoid_monte_carlo.h /$HOME/include
#
gcc -c -I/$HOME/include ellipsoid_monte_carlo.c
if [ $? -ne 0 ]; then
  echo "Errors compiling ellipsoid_monte_carlo.c"
  exit
fi
#
mv ellipsoid_monte_carlo.o ~/libc/$ARCH/ellipsoid_monte_carlo.o
#
echo "Library installed as ~/libc/$ARCH/ellipsoid_monte_carlo.o"
