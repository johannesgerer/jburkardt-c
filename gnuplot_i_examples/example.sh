#!/bin/bash
#
gcc -c -I/$HOME/include example.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling example.c"
  exit
fi
rm compiler.txt
#
gcc example.o ~/libc/$ARCH/gnuplot_i.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading example.o"
  exit
fi
rm example.o
#
mv a.out ~/binc/$ARCH/example
#
echo "Executable installed as ~/binc/$ARCH/example"
