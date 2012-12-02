#!/bin/bash
#
cp condition.h /$HOME/include
#
gcc -c -g -I /$HOME/include condition.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling condition.c"
  exit
fi
rm compiler.txt
#
mv condition.o ~/libc/$ARCH/condition.o
#
echo "Library installed as ~/libc/$ARCH/condition.o"
