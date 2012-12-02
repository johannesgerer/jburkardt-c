#!/bin/bash
#
cp asa103.h /$HOME/include
#
gcc -c -g asa103.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa103.c."
  exit
fi
rm compiler.txt
#
mv asa103.o ~/libc/$ARCH/asa103.o
#
echo "Library installed as ~/libc/$ARCH/asa103.o"
