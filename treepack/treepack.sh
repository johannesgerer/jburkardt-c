#!/bin/bash
#
cp treepack.h /$HOME/include
#
gcc -c -g -I/$HOME/include treepack.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling treepack.c"
  exit
fi
rm compiler.txt
#
mv treepack.o ~/libc/$ARCH/treepack.o
#
echo "Library installed as ~/libc/$ARCH/treepack.o"
