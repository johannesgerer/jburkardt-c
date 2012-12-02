#!/bin/bash
#
cp gm_rule.h /$HOME/include
#
gcc -c -g -I /$HOME/include gm_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gm_rule.c"
  exit
fi
rm compiler.txt
#
mv gm_rule.o ~/libc/$ARCH/gm_rule.o
#
echo "Library installed as ~/libc/$ARCH/gm_rule.o"
