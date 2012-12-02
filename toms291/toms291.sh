#!/bin/bash
#
cp toms291.h /$HOME/include
#
gcc -c -g toms291.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms291.c."
  exit
fi
rm compiler.txt
#
mv toms291.o ~/libc/$ARCH/toms291.o
#
echo "Library installed as ~/libc/$ARCH/toms291.o"
