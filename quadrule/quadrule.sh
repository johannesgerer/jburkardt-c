#!/bin/bash
#
cp quadrule.h /$HOME/include
#
gcc -c -g quadrule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadrule.c."
  exit
fi
rm compiler.txt
#
mv quadrule.o ~/libc/$ARCH/quadrule.o
#
echo "Library installed as ~/libc/$ARCH/quadrule.o"
