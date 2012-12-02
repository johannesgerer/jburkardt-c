#!/bin/bash
#
gcc -c toms178.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms178.c."
  exit
fi
rm compiler.txt
#
gcc toms178.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms178.o -lm"
  exit
fi
#
rm toms178.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/toms178
#
echo "Executable installed as ~/binc/$ARCH/toms178"
