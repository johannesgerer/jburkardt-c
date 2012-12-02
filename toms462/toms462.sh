#!/bin/bash
#
cp toms462.h /$HOME/include
#
gcc -c -g -I /$HOME/include toms462.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms462.c"
  exit
fi
rm compiler.txt
#
mv toms462.o ~/libc/$ARCH/toms462.o
#
echo "Library installed as ~/libc/$ARCH/toms462.o"
