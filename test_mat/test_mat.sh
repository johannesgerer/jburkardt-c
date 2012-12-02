#!/bin/bash
#
cp test_mat.h /$HOME/include
#
gcc -c -g test_mat.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_mat.c."
  exit
fi
rm compiler.txt
#
mv test_mat.o ~/libc/$ARCH/test_mat.o
#
echo "Library installed as ~/libc/$ARCH/test_mat.o"
