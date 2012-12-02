#!/bin/bash
#
cp test_interp_2d.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_interp_2d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_interp_2d.c"
  exit
fi
rm compiler.txt
#
mv test_interp_2d.o ~/libc/$ARCH/test_interp_2d.o
#
echo "Library installed as ~/libc/$ARCH/test_interp_2d.o"
