#!/bin/bash
#
gcc -c legendre_rule_fast.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling legendre_rule_fast.c."
  exit
fi
rm compiler.txt
#
gcc legendre_rule_fast.o -lm
if [ $? -ne 0 ]; then
  echo "Errors loading legendre_rule_fast.c."
  exit
fi
$
rm legendre_rule_fast.o
#
mv a.out ~/binc/$ARCH/legendre_rule_fast
#
echo "Executable installed as ~/binc/$ARCH/legendre_rule_fast"
