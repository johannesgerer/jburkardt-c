#!/bin/bash
#
cp rkf45.h /$HOME/include
#
gcc -c -g rkf45.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rkf45.c."
  exit
fi
rm compiler.txt
#
mv rkf45.o ~/libc/$ARCH/rkf45.o
#
echo "Library installed as ~/libc/$ARCH/rkf45.o"
