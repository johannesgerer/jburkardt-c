#!/bin/bash
#
gcc -c -I/$HOME/include anim.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling anim.c"
  exit
fi
rm compiler.txt
#
gcc anim.o ~/libc/$ARCH/gnuplot_i.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading anim.o"
  exit
fi
rm anim.o
#
mv a.out ~/binc/$ARCH/anim
#
echo "Executable installed as ~/binc/$ARCH/anim"
