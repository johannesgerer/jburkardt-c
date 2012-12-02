#!/bin/bash
#
gcc -c -g diaphony.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling diaphony.c"
  exit
fi
rm compiler.txt
#
gcc diaphony.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading diaphony.o"
  exit
fi
rm diaphony.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/diaphony
#
echo "Executable installed as ~/binc/$ARCH/diaphony"
