#!/bin/bash
#
cp test_optimization.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_optimization.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_optimization.c"
  exit
fi
rm compiler.txt
#
mv test_optimization.o ~/libc/$ARCH/test_optimization.o
#
echo "Library installed as ~/libc/$ARCH/test_optimization.o"
