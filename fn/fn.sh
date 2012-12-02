#!/bin/bash
#
cp fn.h /$HOME/include
#
gcc -c -g -I /$HOME/include fn.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fn.c"
  exit
fi
rm compiler.txt
#
mv fn.o ~/libc/$ARCH/fn.o
#
echo "Library installed as ~/libc/$ARCH/fn.o"
