#!/bin/bash
#
cp test_int_2d.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_int_2d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_int_2d.c"
  exit
fi
rm compiler.txt
#
mv test_int_2d.o ~/libc/$ARCH/test_int_2d.o
#
echo "Library installed as ~/libc/$ARCH/test_int_2d.o"
