#!/bin/bash
#
cp fd1d_bvp.h /$HOME/include
#
gcc -c -g fd1d_bvp.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_bvp.c."
  exit
fi
rm compiler.txt
#
mv fd1d_bvp.o ~/libc/$ARCH/fd1d_bvp.o
#
echo "Library installed as ~/libc/$ARCH/fd1d_bvp.o"
