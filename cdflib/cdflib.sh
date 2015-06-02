#!/bin/bash
#
cp cdflib.h /$HOME/include
#
gcc -c -g -I /$HOME/include cdflib.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cdflib.c"
  exit
fi
rm compiler.txt
#
mv cdflib.o ~/libc/$ARCH/cdflib.o
#
echo "Library installed as ~/libc/$ARCH/cdflib.o"
