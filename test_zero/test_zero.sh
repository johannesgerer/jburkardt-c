#!/bin/bash
#
cp test_zero.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_zero.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_zero.c"
  exit
fi
rm compiler.txt
#
mv test_zero.o ~/libc/$ARCH/test_zero.o
#
echo "Library installed as ~/libc/$ARCH/test_zero.o"
