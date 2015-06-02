#!/bin/bash
#
cp laguerre_test_int.h /$HOME/include
#
gcc -c -g -I /$HOME/include laguerre_test_int.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laguerre_test_int.c"
  exit
fi
rm compiler.txt
#
mv laguerre_test_int.o ~/libc/$ARCH/laguerre_test_int.o
#
echo "Library installed as ~/libc/$ARCH/laguerre_test_int.o"
