#!/bin/bash
#
cp fem1d_pack.h /$HOME/include
#
gcc -c -g -I /$HOME/include fem1d_pack.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_pack.c"
  exit
fi
rm compiler.txt
#
mv fem1d_pack.o ~/libc/$ARCH/fem1d_pack.o
#
echo "Library installed as ~/libc/$ARCH/fem1d_pack.o"
