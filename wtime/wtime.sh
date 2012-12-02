#!/bin/bash
#
cp wtime.h /$HOME/include
#
gcc -c -g wtime.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wtime.c."
  exit
fi
rm compiler.txt
#
mv wtime.o ~/libc/$ARCH/wtime.o
#
echo "Library installed as ~/libc/$ARCH/wtime.o"
