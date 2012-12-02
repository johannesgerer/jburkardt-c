#!/bin/bash
#
cp combo.h /$HOME/include
#
gcc -c -g -I /$HOME/include combo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling combo.c"
  exit
fi
rm compiler.txt
#
mv combo.o ~/libc/$ARCH/combo.o
#
echo "Library installed as ~/libc/$ARCH/combo.o"
