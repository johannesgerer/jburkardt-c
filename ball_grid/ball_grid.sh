#!/bin/bash
#
cp ball_grid.h /$HOME/include
#
gcc -c -g ball_grid.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ball_grid.c."
  exit
fi
rm compiler.txt
#
mv ball_grid.o ~/libc/$ARCH/ball_grid.o
#
echo "Library installed as ~/libc/$ARCH/ball_grid.o"
