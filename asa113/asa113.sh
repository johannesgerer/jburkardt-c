#!/bin/bash
#
cp asa113.h /$HOME/include
#
gcc -c -g asa113.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa113.c."
  exit
fi
rm compiler.txt
#
mv asa113.o ~/libc/$ARCH/asa113.o
#
echo "Library installed as ~/libc/$ARCH/asa113.o"
