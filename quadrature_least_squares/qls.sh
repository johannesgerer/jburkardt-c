#!/bin/bash
#
cp qls.h /$HOME/include
#
gcc -c -I/$HOME/include qls.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qls.c"
  exit
fi
rm compiler.txt
#
mv qls.o ~/libc/$ARCH/qls.o
#
echo "Library installed as ~/libc/$ARCH/qls.o"
