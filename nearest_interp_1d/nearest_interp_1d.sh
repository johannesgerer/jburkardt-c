#!/bin/bash
#
cp nearest_interp_1d.h /$HOME/include
#
gcc -c -g -I /$HOME/include nearest_interp_1d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nearest_interp_1d.c"
  exit
fi
rm compiler.txt
#
mv nearest_interp_1d.o ~/libc/$ARCH/nearest_interp_1d.o
#
echo "Library installed as ~/libc/$ARCH/nearest_interp_1d.o"
