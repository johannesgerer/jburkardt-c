#!/bin/bash
#
cp test_int_hermite.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_int_hermite.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_int_hermite.c"
  exit
fi
rm compiler.txt
#
mv test_int_hermite.o ~/libc/$ARCH/test_int_hermite.o
#
echo "Library installed as ~/libc/$ARCH/test_int_hermite.o"
