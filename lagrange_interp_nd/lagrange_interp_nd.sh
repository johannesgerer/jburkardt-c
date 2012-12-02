#!/bin/bash
#
cp lagrange_interp_nd.h /$HOME/include
#
gcc -c -g -I/$HOME/include lagrange_interp_nd.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lagrange_interp_nd.c"
  exit
fi
rm compiler.txt
#
mv lagrange_interp_nd.o ~/libc/$ARCH/lagrange_interp_nd.o
#
echo "Library installed as ~/libc/$ARCH/lagrange_interp_nd.o"
