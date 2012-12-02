#!/bin/bash
#
cp asa109.h /$HOME/include
#
gcc -c -g asa109.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa109.c."
  exit
fi
rm compiler.txt
#
mv asa109.o ~/libc/$ARCH/asa109.o
#
echo "Library installed as ~/libc/$ARCH/asa109.o"
