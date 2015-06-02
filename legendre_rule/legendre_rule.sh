#!/bin/bash
#
gcc -c -g legendre_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling legendre_rule.c"
  exit
fi
rm compiler.txt
#
gcc legendre_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading legendre_rule.o"
  exit
fi
rm legendre_rule.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/legendre_rule
#
echo "Executable installed as ~/binc/$ARCH/legendre_rule"
