#!/bin/bash
#
cp rbf_interp_2d.h /$HOME/include
#
gcc -c -g -I/$HOME/include rbf_interp_2d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rbf_interp_2d.c"
  exit
fi
rm compiler.txt
#
mv rbf_interp_2d.o ~/libc/$ARCH/rbf_interp_2d.o
#
echo "Library installed as ~/libc/$ARCH/rbf_interp_2d.o"
