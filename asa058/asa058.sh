#!/bin/bash
#
cp asa058.h /$HOME/include
#
gcc -c -g asa058.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa058.c."
  exit
fi
rm compiler.txt
#
mv asa058.o ~/libc/$ARCH/asa058.o
#
echo "Library installed as ~/libc/$ARCH/asa058.o"
