#!/bin/bash
#
cp test_int.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_int.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_int.c."
  exit
fi
rm compiler.txt
#
mv test_int.o ~/libc/$ARCH/test_int.o
#
echo "Library installed as ~/libc/$ARCH/test_int.o"
