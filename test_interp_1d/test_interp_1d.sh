#!/bin/bash
#
cp test_interp_1d.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_interp_1d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_interp_1d.c"
  exit
fi
rm compiler.txt
#
mv test_interp_1d.o ~/libc/$ARCH/test_interp_1d.o
#
echo "Library installed as ~/libc/$ARCH/test_interp_1d.o"
