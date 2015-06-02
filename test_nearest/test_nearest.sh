#!/bin/bash
#
gcc -c -g test_nearest.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_nearest.c"
  exit
fi
rm compiler.txt
#
gcc test_nearest.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_nearest.o"
  exit
fi
rm test_nearest.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/test_nearest
#
echo "Executable installed as ~/binc/$ARCH/test_nearest"
