#!/bin/bash
#
cp fem2d_pack.h /$HOME/include
#
gcc -c -g -I/$HOME/include fem2d_pack.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_pack.c"
  exit
fi
rm compiler.txt
#
mv fem2d_pack.o ~/libc/$ARCH/fem2d_pack.o
#
echo "Library installed as ~/libc/$ARCH/fem2d_pack.o"
