#!/bin/bash
#
gcc -c -g fd1d_burgers_leap.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling fd1d_burgers_leap.c"
  exit
fi
rm compiler.txt
#
gcc fd1d_burgers_leap.o
if [ $? -ne 0 ]; then
  echo "Errors while loading fd1d_burgers_leap.o"
  exit
fi
rm fd1d_burgers_leap.o
#
mv a.out ~/binc/$ARCH/fd1d_burgers_leap
#
echo "Program installed as ~/binc/$ARCH/fd1d_burgers_leap"
