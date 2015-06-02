#!/bin/bash
#
gcc -c -g monte_carlo_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling monte_carlo_rule.c"
  exit
fi
rm compiler.txt
#
gcc monte_carlo_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading monte_carlo_rule.o"
  exit
fi
rm monte_carlo_rule.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/monte_carlo_rule
#
echo "Executable installed as ~/binc/$ARCH/monte_carlo_rule"
