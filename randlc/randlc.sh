#!/bin/bash
#
cp randlc.h /$HOME/include
#
gcc -c -g randlc.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling randlc.c."
  exit
fi
rm compiler.txt
#
mv randlc.o ~/libc/$ARCH/randlc.o
#
echo "Library installed as ~/libc/$ARCH/randlc.o"
