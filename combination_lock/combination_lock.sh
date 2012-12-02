#!/bin/bash
#
cp combination_lock.h /$HOME/include
#
gcc -c -g -I /$HOME/include combination_lock.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling combination_lock.c"
  exit
fi
rm compiler.txt
#
mv combination_lock.o ~/libc/$ARCH/combination_lock.o
#
echo "Library installed as ~/libc/$ARCH/combination_lock.o"
