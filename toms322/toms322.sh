#!/bin/bash
#
cp toms322.h /$HOME/include
#
gcc -c -g toms322.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms322.c."
  exit
fi
rm compiler.txt
#
mv toms322.o ~/libc/$ARCH/toms322.o
#
echo "Library installed as ~/libc/$ARCH/toms322.o"
