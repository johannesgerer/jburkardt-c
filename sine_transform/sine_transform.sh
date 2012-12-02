#!/bin/bash
#
cp sine_transform.h /$HOME/include
#
gcc -c -g -I /$HOME/include sine_transform.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sine_transform.c"
  exit
fi
rm compiler.txt
#
mv sine_transform.o ~/libc/$ARCH/sine_transform.o
#
echo "Library installed as ~/libc/$ARCH/sine_transform.o"
