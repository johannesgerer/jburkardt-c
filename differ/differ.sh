#!/bin/bash
#
cp differ.h /$HOME/include
#
gcc -c -g -I/$HOME/include differ.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling differ.c"
  exit
fi
rm compiler.txt
#
mv differ.o ~/libc/$ARCH/differ.o
#
echo "Library installed as ~/libc/$ARCH/differ.o"
