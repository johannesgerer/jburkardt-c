#!/bin/bash
#
cp asa183.h /$HOME/include
#
gcc -c -g asa183.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa183.c."
  exit
fi
rm compiler.txt
#
mv asa183.o ~/libc/$ARCH/asa183.o
#
echo "Library installed as ~/libc/$ARCH/asa183.o"
