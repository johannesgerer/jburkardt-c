#!/bin/bash
#
cp barycentric_interp_1d.h /$HOME/include
#
gcc -c -g -I/$HOME/include barycentric_interp_1d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling barycentric_interp_1d.c"
  exit
fi
rm compiler.txt
#
mv barycentric_interp_1d.o ~/libc/$ARCH/barycentric_interp_1d.o
#
echo "Library installed as ~/libc/$ARCH/barycentric_interp_1d.o"
