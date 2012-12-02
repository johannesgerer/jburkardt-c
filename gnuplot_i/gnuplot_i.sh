#!/bin/bash
#
cp gnuplot_i.h /$HOME/include
cp gnuplot_i.hpp /$HOME/include
#
gcc -c gnuplot_i.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gnuplot_i.c"
  exit
fi
rm compiler.txt
#
mv gnuplot_i.o ~/libc/$ARCH/gnuplot_i.o
#
echo "Library installed as ~/libc/$ARCH/gnuplot_i.o"
