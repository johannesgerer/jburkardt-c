#!/bin/bash
#
gcc -c -g complex_numbers.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling complex_numbers.c"
  exit
fi
rm compiler.txt
#
gcc complex_numbers.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading complex_numbers.o"
  exit
fi
rm complex_numbers.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/complex_numbers
#
echo "Executable installed as ~/binc/$ARCH/complex_numbers"
