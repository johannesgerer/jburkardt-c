#!/bin/bash
#
cp asa136.h /$HOME/include
#
gcc -c -g asa136.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa136.c."
  exit
fi
rm compiler.txt
#
mv asa136.o ~/libc/$ARCH/asa136.o
#
echo "Library installed as ~/libc/$ARCH/asa136.o"
