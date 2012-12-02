#!/bin/bash
#
cp asa152.h /$HOME/include
#
gcc -c -g asa152.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa152.c."
  exit
fi
rm compiler.txt
#
mv asa152.o ~/libc/$ARCH/asa152.o
#
echo "Library installed as ~/libc/$ARCH/asa152.o"
