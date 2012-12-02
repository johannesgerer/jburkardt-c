#!/bin/bash
#
gcc -c -g -I/$HOME/include ply_to_obj.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ply_to_obj.c."
  exit
fi
rm compiler.txt
#
gcc ply_to_obj.o $HOME/libc/$ARCH/ply_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ply_to_obj.o + ply_io.o."
  exit
fi
#
rm ply_to_obj.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/ply_to_obj
#
echo "Executable installed as ~/binc/$ARCH/ply_to_obj"
