#!/bin/bash
#
gcc -c -I/$HOME/include obj_to_ply.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling obj_to_ply.c."
  exit
fi
rm compiler.txt
#
gcc obj_to_ply.o $HOME/libc/$ARCH/ply_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading obj_to_ply.o + ply_io.o."
  exit
fi
#
rm obj_to_ply.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/obj_to_ply
#
echo "Program installed as ~/binc/$ARCH/obj_to_ply"
