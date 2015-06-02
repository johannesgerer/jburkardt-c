#!/bin/bash
#
cp backtrack_binary_rc.h /$HOME/include
#
gcc -c -g -I/$HOME/include backtrack_binary_rc.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling backtrack_binary_rc.c"
  exit
fi
rm compiler.txt
#
mv backtrack_binary_rc.o ~/libc/$ARCH/backtrack_binary_rc.o
#
echo "Library installed as ~/libc/$ARCH/backtrack_binary_rc.o"
