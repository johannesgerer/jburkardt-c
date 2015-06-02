#!/bin/bash
#
gcc -c -g truncated_normal_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling truncated_normal_rule.c"
  exit
fi
#
gcc truncated_normal_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors loading truncated_normal_rule.c"
  exit
fi
#
rm truncated_normal_rule.o
mv a.out ~/binc/$ARCH/truncated_normal_rule
#
echo "Executable installed as ~/binc/$ARCH/truncated_normal_rule."
