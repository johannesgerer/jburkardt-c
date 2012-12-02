#!/bin/bash
#
cp toms179.h /$HOME/include
#
gcc -c -g toms179.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms179.c."
  exit
fi
rm compiler.txt
#
mv toms179.o ~/libc/$ARCH/toms179.o
#
echo "Library installed as ~/libc/$ARCH/toms179.o"
