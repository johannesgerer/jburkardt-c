#!/bin/bash
#
cp asa241.h /$HOME/include
#
gcc -c -g asa241.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa241.c."
  exit
fi
rm compiler.txt
#
mv asa241.o ~/libc/$ARCH/asa241.o
#
echo "Library installed as ~/libc/$ARCH/asa241.o"
