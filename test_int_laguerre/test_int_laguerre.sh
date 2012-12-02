#!/bin/bash
#
cp test_int_laguerre.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_int_laguerre.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_int_laguerre.c"
  exit
fi
rm compiler.txt
#
mv test_int_laguerre.o ~/libc/$ARCH/test_int_laguerre.o
#
echo "Library installed as ~/libc/$ARCH/test_int_laguerre.o"
