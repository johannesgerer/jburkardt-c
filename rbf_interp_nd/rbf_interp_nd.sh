#!/bin/bash
#
cp rbf_interp_nd.h /$HOME/include
#
gcc -c -g -I /$HOME/include rbf_interp_nd.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rbf_interp_nd.c"
  exit
fi
rm compiler.txt
#
mv rbf_interp_nd.o ~/libc/$ARCH/rbf_interp_nd.o
#
echo "Library installed as ~/libc/$ARCH/rbf_interp_nd.o"
