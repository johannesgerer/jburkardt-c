#!/bin/bash
#
cp asa147.h /$HOME/include
#
gcc -c -g asa147.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa147.c."
  exit
fi
rm compiler.txt
#
mv asa147.o ~/libc/$ARCH/asa147.o
#
echo "Library installed as ~/libc/$ARCH/asa147.o"
