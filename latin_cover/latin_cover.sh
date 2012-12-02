#!/bin/bash
#
cp latin_cover.h /$HOME/include
#
gcc -c -g latin_cover.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_cover.c"
  exit
fi
rm compiler.txt
#
mv latin_cover.o ~/libc/$ARCH/latin_cover.o
#
echo "Library installed as ~/libc/$ARCH/latin_cover.o"
