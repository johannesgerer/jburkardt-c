#!/bin/bash
#
cp qwgw.h /$HOME/include
#
gcc -c -g -I/$HOME/include qwgw.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qwgw.c"
  exit
fi
rm compiler.txt
#
mv qwgw.o ~/libc/$ARCH/qwgw.o
#
echo "Library installed as ~/libc/$ARCH/qwgw.o"
