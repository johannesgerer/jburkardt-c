#!/bin/bash
#
cp asa006.h /$HOME/include
#
gcc -c -g asa006.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa006.c."
  exit
fi
rm compiler.txt
#
mv asa006.o ~/libc/$ARCH/asa006.o
#
echo "Library installed as ~/libc/$ARCH/asa006.o"
