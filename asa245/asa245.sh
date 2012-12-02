#!/bin/bash
#
cp asa245.h /$HOME/include
#
gcc -c -g asa245.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa245.c."
  exit
fi
rm compiler.txt
#
mv asa245.o ~/libc/$ARCH/asa245.o
#
echo "Library installed as ~/libc/$ARCH/asa245.o"
