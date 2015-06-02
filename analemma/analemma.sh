#!/bin/bash
#
gcc -c -g analemma.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling analemma.c"
  exit
fi
rm compiler.txt
#
gcc analemma.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading analemma.o"
  exit
fi
#
rm analemma.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/analemma
#
echo "Executable installed as ~/binc/$ARCH/analemma"
