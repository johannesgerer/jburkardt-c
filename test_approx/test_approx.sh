#!/bin/bash
#
cp test_approx.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_approx.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_approx.c"
  exit
fi
rm compiler.txt
#
mv test_approx.o ~/libc/$ARCH/test_approx.o
#
echo "Library installed as ~/libc/$ARCH/test_approx.o"
