#!/bin/bash
#
gcc -c -g -I/$HOME/include ply_to_iv.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ply_to_iv.c."
  exit
fi
rm compiler.txt
#
gcc ply_to_iv.o $HOME/libc/$ARCH/ply_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ply_to_iv.o + ply_io.o."
  exit
fi
#
rm ply_to_iv.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/ply_to_iv
#
echo "Executable installed as ~/binc/$ARCH/ply_to_iv"
