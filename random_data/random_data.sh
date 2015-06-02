#!/bin/bash
#
cp random_data.h /$HOME/include
#
gcc -c -g random_data.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling random_data.c"
  exit
fi
rm compiler.txt
#
mv random_data.o ~/libc/$ARCH/random_data.o
#
echo "Library installed as ~/libc/$ARCH/random_data.o"
