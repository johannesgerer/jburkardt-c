#!/bin/bash
#
cp asa299.h /$HOME/include
#
gcc -c -g asa299.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa299.c."
  exit
fi
rm compiler.txt
#
mv asa299.o ~/libc/$ARCH/asa299.o
#
echo "Library installed as ~/libc/$ARCH/asa299.o"
