#!/bin/bash
#
cp rbf_interp_1d.h /$HOME/include
#
gcc -c -g -I/$HOME/include rbf_interp_1d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rbf_interp_1d.c"
  exit
fi
rm compiler.txt
#
mv rbf_interp_1d.o ~/libc/$ARCH/rbf_interp_1d.o
#
echo "Library installed as ~/libc/$ARCH/rbf_interp_1d.o"
