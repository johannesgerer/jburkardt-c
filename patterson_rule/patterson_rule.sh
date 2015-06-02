#!/bin/bash
#
gcc -c patterson_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling patterson_rule.c"
  exit
fi
rm compiler.txt
#
gcc patterson_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading patterson_rule.o"
  exit
fi
rm patterson_rule.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/patterson_rule
#
echo "Executable installed as ~/binc/$ARCH/patterson_rule"
