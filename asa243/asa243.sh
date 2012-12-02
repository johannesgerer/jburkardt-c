#!/bin/bash
#
cp asa243.h /$HOME/include
#
gcc -c -g asa243.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa243.c."
  exit
fi
rm compiler.txt
#
mv asa243.o ~/libc/$ARCH/asa243.o
#
echo "Library installed as ~/libc/$ARCH/asa243.o"
