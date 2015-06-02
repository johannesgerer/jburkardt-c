#!/bin/bash
#
cp circle_rule.h /$HOME/include
#
gcc -c -g -I /$HOME/include circle_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_rule.c"
  exit
fi
rm compiler.txt
#
mv circle_rule.o ~/libc/$ARCH/circle_rule.o
#
echo "Library installed as ~/libc/$ARCH/circle_rule.o"
