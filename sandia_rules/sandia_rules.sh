#!/bin/bash
#
cp sandia_rules.h /$HOME/include
#
gcc -c -g sandia_rules.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_rules.c."
  exit
fi
rm compiler.txt
#
mv sandia_rules.o ~/libc/$ARCH/sandia_rules.o
#
echo "Library installed as ~/libc/$ARCH/sandia_rules.o"
