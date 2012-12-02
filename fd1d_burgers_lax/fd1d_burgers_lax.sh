#!/bin/bash
#
gcc -c -g fd1d_burgers_lax.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling fd1d_burgers_lax.c"
  exit
fi
rm compiler.txt
#
gcc fd1d_burgers_lax.o
if [ $? -ne 0 ]; then
  echo "Errors while loading fd1d_burgers_lax.o"
  exit
fi
rm fd1d_burgers_lax.o
#
mv a.out ~/binc/$ARCH/fd1d_burgers_lax
#
echo "Program installed as ~/binc/$ARCH/fd1d_burgers_lax"
