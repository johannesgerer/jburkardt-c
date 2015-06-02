#!/bin/bash
#
cp toms886.h /$HOME/include
#
gcc -c -g -I/$HOME/include toms886.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms886.c"
  exit
fi
rm compiler.txt
#
mv toms886.o ~/libc/$ARCH/toms886.o
#
echo "Library installed as ~/libc/$ARCH/toms886.o"
