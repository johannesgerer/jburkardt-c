#!/bin/bash
#
cp asa053.h /$HOME/include
#
gcc -c -I/$HOME/include asa053.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa053.c"
  exit
fi
rm compiler.txt
#
mv asa053.o ~/libc/$ARCH/asa053.o
#
echo "Library installed as ~/libc/$ARCH/asa053.o"
