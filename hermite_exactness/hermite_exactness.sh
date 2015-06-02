#!/bin/bash
#
gcc -c -g hermite_exactness.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_exactness.c"
  exit
fi
rm compiler.txt
#
gcc hermite_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hermite_exactness.o"
  exit
fi
rm hermite_exactness.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/hermite_exactness
#
echo "Executable installed as ~/binc/$ARCH/hermite_exactness"
