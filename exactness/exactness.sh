#!/bin/bash
#
cp exactness.h /$HOME/include
#
gcc -c -I/$HOME/include exactness.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling exactness.c"
  exit
fi
rm compiler.txt
#
mv exactness.o ~/libc/$ARCH/exactness.o
#
echo "Library installed as ~/libc/$ARCH/exactness.o"
