#!/bin/bash
#
gcc -c pyramid_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_rule.c"
  exit
fi
#
gcc pyramid_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pyramid_rule.o"
  exit
fi
rm pyramid_rule.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/pyramid_rule
#
echo "Executable installed as ~/binc/$ARCH/pyramid_rule"
