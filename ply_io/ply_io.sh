#!/bin/bash
#
cp ply_io.h /$HOME/include
#
gcc -c -g -I/$HOME/include ply_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ply_io.c."
  exit
fi
rm compiler.txt
#
mv ply_io.o ~/libc/$ARCH/ply_io.o
#
echo "Library installed as ~/libc/$ARCH/ply_io.o"
