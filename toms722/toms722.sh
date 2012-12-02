#!/bin/bash
#
cp toms722.h /$HOME/include
#
gcc -c -g toms722.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms722.c."
  exit
fi
rm compiler.txt
#
mv toms722.o ~/libc/$ARCH/toms722.o
#
echo "Library installed as ~/libc/$ARCH/toms722.o"
