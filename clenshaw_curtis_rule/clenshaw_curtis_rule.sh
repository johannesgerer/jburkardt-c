#!/bin/bash
#
gcc -c -g clenshaw_curtis_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling clenshaw_curtis_rule.c"
  exit
fi
rm compiler.txt
#
gcc clenshaw_curtis_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading clenshaw_curtis_rule.o"
  exit
fi
rm clenshaw_curtis_rule.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/clenshaw_curtis_rule
#
echo "Executable installed as ~/binc/$ARCH/clenshaw_curtis_rule"
