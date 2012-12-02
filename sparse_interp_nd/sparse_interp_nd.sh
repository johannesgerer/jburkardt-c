#!/bin/bash
#
cp sparse_interp_nd.h /$HOME/include
#
gcc -c -g -I/$HOME/include sparse_interp_nd.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_interp_nd.c"
  exit
fi
rm compiler.txt
#
mv sparse_interp_nd.o ~/libc/$ARCH/sparse_interp_nd.o
#
echo "Library installed as ~/libc/$ARCH/sparse_interp_nd.o"
