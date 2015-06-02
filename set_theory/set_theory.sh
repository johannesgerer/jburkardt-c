#!/bin/bash
#
cp set_theory.h /$HOME/include
#
gcc -c -g -I /$HOME/include set_theory.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling set_theory.c"
  exit
fi
rm compiler.txt
#
mv set_theory.o ~/libc/$ARCH/set_theory.o
#
echo "Library installed as ~/libc/$ARCH/set_theory.o"
