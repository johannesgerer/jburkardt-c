#!/bin/bash
#
gcc -c -I/$HOME/include plato_ply.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling plato_ply.c."
  exit
fi
rm compiler.txt
#
gcc plato_ply.o $HOME/libc/$ARCH/ply_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading plato_ply.o + ply_io.o."
  exit
fi
#
rm plato_ply.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/plato_ply
#
echo "Executable installed as ~/binc/$ARCH/plato_ply"
