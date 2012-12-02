#!/bin/bash
#
cp test_values.h /$HOME/include
#
gcc -c -g test_values.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_values.c."
  exit
fi
rm compiler.txt
#
mv test_values.o ~/libc/$ARCH/test_values.o
#
echo "Library installed as ~/libc/$ARCH/test_values.o"
