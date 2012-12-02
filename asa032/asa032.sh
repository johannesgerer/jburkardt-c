#!/bin/bash
#
cp asa032.h /$HOME/include
#
gcc -c -g asa032.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa032.c."
  exit
fi
rm compiler.txt
#
mv asa032.o ~/libc/$ARCH/asa032.o
#
echo "Library installed as ~/libc/$ARCH/asa032.o"
