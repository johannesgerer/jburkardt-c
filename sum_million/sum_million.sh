#!/bin/bash
#
gcc -c sum_million.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sum_million.c"
  exit
fi
rm compiler.txt
#
gcc sum_million.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sum_million.o"
  exit
fi
rm sum_million.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/sum_million
#
echo "Program installed as ~/binc/$ARCH/sum_million"
