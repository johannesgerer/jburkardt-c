#!/bin/bash
#
gcc -c damped_sine.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling damped_sine.c"
  exit
fi
rm compiler.txt
#
gcc damped_sine.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking damped_sine.o"
  exit
fi
#
rm damped_sine.o
mv a.out ~/binc/$ARCH/damped_sine
#
echo "Executable installed as ~/binc/$ARCH/damped_sine"
