#!/bin/bash
#
cp haar.h /$HOME/include
#
gcc -c -g -I /$HOME/include haar.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling haar.c."
  exit
fi
rm compiler.txt
#
mv haar.o ~/libc/$ARCH/haar.o
#
echo "Library installed as ~/libc/$ARCH/haar.o"
