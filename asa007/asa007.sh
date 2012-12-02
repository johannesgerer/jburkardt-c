#!/bin/bash
#
cp asa007.h /$HOME/include
#
gcc -c -g asa007.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa007.c."
  exit
fi
rm compiler.txt
#
mv asa007.o ~/libc/$ARCH/asa007.o
#
echo "Library installed as ~/libc/$ARCH/asa007.o"
