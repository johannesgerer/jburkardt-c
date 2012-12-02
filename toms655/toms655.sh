#!/bin/bash
#
cp toms655.h /$HOME/include
#
gcc -c -g toms655.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms655.c."
  exit
fi
rm compiler.txt
#
mv toms655.o ~/libc/$ARCH/toms655.o
#
echo "Library installed as ~/libc/$ARCH/toms655.o"
