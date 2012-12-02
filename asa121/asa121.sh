#!/bin/bash
#
cp asa121.h /$HOME/include
#
gcc -c -g asa121.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa121.c."
  exit
fi
rm compiler.txt
#
mv asa121.o ~/libc/$ARCH/asa121.o
#
echo "Library installed as ~/libc/$ARCH/asa121.o"
