#!/bin/bash
#
cp wishart.h /$HOME/include
#
gcc -c -g -I/$HOME/include wishart.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wishart.c"
  exit
fi
rm compiler.txt
#
mv wishart.o ~/libc/$ARCH/wishart.o
#
echo "Library installed as ~/libc/$ARCH/wishart.o"
