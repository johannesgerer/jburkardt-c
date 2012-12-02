#!/bin/bash
#
cp bernstein.h /$HOME/include
#
gcc -c -g -I /$HOME/include bernstein.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bernstein.c"
  exit
fi
rm compiler.txt
#
mv bernstein.o ~/libc/$ARCH/bernstein.o
#
echo "Library installed as ~/libc/$ARCH/bernstein.o"
