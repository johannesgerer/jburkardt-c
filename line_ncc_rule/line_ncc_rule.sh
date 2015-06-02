#!/bin/bash
#
cp line_ncc_rule.h /$HOME/include
#
gcc -c -g -I/$HOME/include line_ncc_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling line_ncc_rule.c"
  exit
fi
rm compiler.txt
#
mv line_ncc_rule.o ~/libc/$ARCH/line_ncc_rule.o
#
echo "Library installed as ~/libc/$ARCH/line_ncc_rule.o"
