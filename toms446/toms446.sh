#!/bin/bash
#
cp toms446.h /$HOME/include
#
gcc -c -g -I /$HOME/include toms446.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms446.c"
  exit
fi
rm compiler.txt
#
mv toms446.o ~/libc/$ARCH/toms446.o
#
echo "Library installed as ~/libc/$ARCH/toms446.o"
