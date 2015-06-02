#!/bin/bash
#
gcc -c -g hermite_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_rule.c"
  exit
fi
rm compiler.txt
#
gcc hermite_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hermite_rule.o"
  exit
fi
rm hermite_rule.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/hermite_rule
#
echo "Executable installed as ~/binc/$ARCH/hermite_rule"
