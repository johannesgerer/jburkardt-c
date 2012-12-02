#!/bin/bash
#
cp linpack_d.h /$HOME/include
#
gcc -c -g -I/$HOME/include linpack_d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_d.c."
  exit
fi
rm compiler.txt
#
mv linpack_d.o ~/libc/$ARCH/linpack_d.o
#
echo "Library installed as ~/libc/$ARCH/linpack_d.o"
