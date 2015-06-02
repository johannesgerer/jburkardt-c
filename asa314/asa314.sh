#!/bin/bash
#
cp asa314.h /$HOME/include
#
gcc -c -g -I/$HOME/include asa314.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa314.c"
  exit
fi
rm compiler.txt
#
mv asa314.o ~/libc/$ARCH/asa314.o
#
echo "Library installed as ~/libc/$ARCH/asa314.o"
