#!/bin/bash
#
cp owens.h /$HOME/include
#
gcc -c -g owens.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling owens.c."
  exit
fi
rm compiler.txt
#
mv owens.o ~/libc/$ARCH/owens.o
#
echo "Library installed as ~/libc/$ARCH/owens.o"
