#!/bin/bash
#
cp asa091.h /$HOME/include
#
gcc -c -g asa091.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa091.c."
  exit
fi
rm compiler.txt
#
mv asa091.o ~/libc/$ARCH/asa091.o
#
echo "Library installed as ~/libc/$ARCH/asa091.o"
