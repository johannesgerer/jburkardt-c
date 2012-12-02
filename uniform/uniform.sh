#!/bin/bash
#
cp uniform.h /$HOME/include
#
gcc -c -g uniform.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling uniform.c."
  exit
fi
rm compiler.txt
#
mv uniform.o ~/libc/$ARCH/uniform.o
#
echo "Library installed as ~/libc/$ARCH/uniform.o"
