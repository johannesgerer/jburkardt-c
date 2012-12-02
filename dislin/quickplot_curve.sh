#!/bin/bash
#
gcc -c -I/usr/local/dislin quickplot_curve.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quickplot_curve.c"
  exit
fi
rm compiler.txt
#
gcc quickplot_curve.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quickplot_curve.o."
  exit
fi
#
rm quickplot_curve.o
#
mv a.out quickplot_curve
./quickplot_curve > quickplot_curve_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quickplot_curve."
  exit
fi
rm quickplot_curve
#
echo "Program output written to quickplot_curve_output.txt"
