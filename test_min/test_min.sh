#!/bin/bash
#
cp test_min.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_min.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_min.c"
  exit
fi
rm compiler.txt
#
mv test_min.o ~/libc/$ARCH/test_min.o
#
echo "Library installed as ~/libc/$ARCH/test_min.o"
