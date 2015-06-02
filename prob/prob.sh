#!/bin/bash
#
cp prob.h /$HOME/include
#
gcc -c -g -I/$HOME/include prob.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling prob.c"
  exit
fi
rm compiler.txt
#
mv prob.o ~/libc/$ARCH/prob.o
#
echo "Library installed as ~/libc/$ARCH/prob.o"
