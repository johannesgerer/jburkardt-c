#!/bin/bash
#
cp test_interp.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_interp.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_interp.c"
  exit
fi
rm compiler.txt
#
mv test_interp.o ~/libc/$ARCH/test_interp.o
#
echo "Library installed as ~/libc/$ARCH/test_interp.o"
