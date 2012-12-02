#!/bin/bash
#
cp csparse.h /$HOME/include
#
gcc -c -g csparse.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling csparse.c."
  exit
fi
rm compiler.txt
#
mv csparse.o ~/libc/$ARCH/csparse.o
#
echo "Library installed as ~/libc/$ARCH/csparse.o"
