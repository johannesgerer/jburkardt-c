#!/bin/bash
#
cp test_interp_nd.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_interp_nd.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_interp_nd.c"
  exit
fi
rm compiler.txt
#
mv test_interp_nd.o ~/libc/$ARCH/test_interp_nd.o
#
echo "Library installed as ~/libc/$ARCH/test_interp_nd.o"
