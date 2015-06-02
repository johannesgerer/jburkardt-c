#!/bin/bash
#
cp toms097.h /$HOME/include
#
gcc -c -g -I/$HOME/include toms097.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms097.c"
  exit
fi
rm compiler.txt
#
mv toms097.o ~/libc/$ARCH/toms097.o
#
echo "Library installed as ~/libc/$ARCH/toms097.o"
