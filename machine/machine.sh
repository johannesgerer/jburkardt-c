#!/bin/bash
#
cp machine.h /$HOME/include
#
gcc -c -g machine.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling machine.c."
  exit
fi
rm compiler.txt
#
mv machine.o ~/libc/$ARCH/machine.o
#
echo "Library installed as ~/libc/$ARCH/machine.o"
