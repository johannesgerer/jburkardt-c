#!/bin/bash
#
cp asa111.h /$HOME/include
#
gcc -c -g asa111.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa111.c."
  exit
fi
rm compiler.txt
#
mv asa111.o ~/libc/$ARCH/asa111.o
#
echo "Library installed as ~/libc/$ARCH/asa111.o"
