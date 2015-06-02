#!/bin/bash
#
gcc -c tsp_brute.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tsp_brute.c"
  exit
fi
rm compiler.txt
#
gcc tsp_brute.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tsp_brute.o"
  exit
fi
rm tsp_brute.o
#
mv a.out ~/binc/$ARCH/tsp_brute
#
echo "Executable installed as ~/binc/$ARCH/tsp_brute"
