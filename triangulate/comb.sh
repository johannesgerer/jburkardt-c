#!/bin/bash
#
gcc -c comb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling comb.c."
  exit
fi
rm compiler.txt
#
gcc comb.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading comb.o"
  exit
fi
#
rm comb.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/comb
#
echo "Executable installed as ~/binc/$ARCH/comb"
