#!/bin/bash
#
gcc -c laguerre_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laguerre_rule.c"
  exit
fi
rm compiler.txt
#
gcc laguerre_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading laguerre_rule.o"
  exit
fi
rm laguerre_rule.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/laguerre_rule
#
echo "Executable installed as ~/binc/$ARCH/laguerre_rule"
