#!/bin/bash
#
cp index.h /$HOME/include
#
gcc -c -g -I/$HOME/include index.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling index.c"
  exit
fi
rm compiler.txt
#
mv index.o ~/libc/$ARCH/index.o
#
echo "Library installed as ~/libc/$ARCH/index.o"
