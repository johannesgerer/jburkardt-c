#!/bin/bash
#
cp hyperball_integrals.h /$HOME/include
#
gcc -c -g -I/$HOME/include hyperball_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hyperball_integrals.c"
  exit
fi
rm compiler.txt
#
mv hyperball_integrals.o ~/libc/$ARCH/hyperball_integrals.o
#
echo "Library installed as ~/libc/$ARCH/hyperball_integrals.o"
