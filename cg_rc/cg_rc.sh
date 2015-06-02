#!/bin/bash
#
cp cg_rc.h /$HOME/include
#
gcc -c -g -I/$HOME/include cg_rc.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cg_rc.c"
  exit
fi
rm compiler.txt
#
mv cg_rc.o ~/libc/$ARCH/cg_rc.o
#
echo "Library installed as ~/libc/$ARCH/cg_rc.o"
