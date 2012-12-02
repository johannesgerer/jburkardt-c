#!/bin/bash
#
gcc -c -g -I $HOME/include memory_test.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling memory_test.c."
  exit
fi
rm compiler.txt
#
gcc memory_test.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading memory_test.o."
  exit
fi
#
rm memory_test.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/memory_test
#
echo "Executable installed as ~/binc/$ARCH/memory_test"
