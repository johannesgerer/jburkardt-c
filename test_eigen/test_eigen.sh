#!/bin/bash
#
cp test_eigen.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_eigen.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_eigen.c"
  exit
fi
rm compiler.txt
#
mv test_eigen.o ~/libc/$ARCH/test_eigen.o
#
echo "Library installed as ~/libc/$ARCH/test_eigen.o"
