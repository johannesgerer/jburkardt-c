#!/bin/bash
#
gcc -c -g -I$HOME/include power_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling power_rule.c"
  exit
fi
rm compiler.txt
#
gcc power_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading power_rule.o."
  exit
fi
#
rm power_rule.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/power_rule
#
echo "Executable installed as ~/binc/$ARCH/power_rule"
