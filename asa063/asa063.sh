#!/bin/bash
#
cp asa063.h /$HOME/include
#
gcc -c -g asa063.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa063.c."
  exit
fi
rm compiler.txt
#
mv asa063.o ~/libc/$ARCH/asa063.o
#
echo "Library installed as ~/libc/$ARCH/asa063.o"
