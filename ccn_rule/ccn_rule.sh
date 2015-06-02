#!/bin/bash
#
gcc -c ccn_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ccn_rule.c"
  exit
fi
rm compiler.txt
#
gcc ccn_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ccn_rule.o"
  exit
fi
rm ccn_rule.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/ccn_rule
#
echo "Executable installed as ~/binc/$ARCH/ccn_rule"
