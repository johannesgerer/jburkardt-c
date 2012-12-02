#!/bin/bash
#
cp asa310.h /$HOME/include
#
gcc -c -g asa310.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa310.c."
  exit
fi
rm compiler.txt
#
mv asa310.o ~/libc/$ARCH/asa310.o
#
echo "Library installed as ~/libc/$ARCH/asa310.o"
