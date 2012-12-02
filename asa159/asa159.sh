#!/bin/bash
#
cp asa159.h /$HOME/include
#
gcc -c -g asa159.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa159.c."
  exit
fi
rm compiler.txt
#
mv asa159.o ~/libc/$ARCH/asa159.o
#
echo "Library installed as ~/libc/$ARCH/asa159.o"
