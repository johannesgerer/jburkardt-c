#!/bin/bash
#
cp llsq.h /$HOME/include
#
gcc -c -g -I /$HOME/include llsq.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling llsq.c"
  exit
fi
rm compiler.txt
#
mv llsq.o ~/libc/$ARCH/llsq.o
#
echo "Library installed as ~/libc/$ARCH/llsq.o"
