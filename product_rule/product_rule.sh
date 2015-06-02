#!/bin/bash
#
gcc -c -g -I$HOME/include product_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling product_rule.c"
  exit
fi
rm compiler.txt
#
gcc product_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading product_rule.o"
  exit
fi
rm product_rule.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/product_rule
#
echo "Executable installed as ~/binc/$ARCH/product_rule"
