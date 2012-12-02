#!/bin/bash
#
cp asa144.h /$HOME/include
#
gcc -c -g asa144.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa144.c."
  exit
fi
rm compiler.txt
#
mv asa144.o ~/libc/$ARCH/asa144.o
#
echo "Library installed as ~/libc/$ARCH/asa144.o"
