#!/bin/bash
#
cp asa239.h /$HOME/include
#
gcc -c -g asa239.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa239.c."
  exit
fi
rm compiler.txt
#
mv asa239.o ~/libc/$ARCH/asa239.o
#
echo "Library installed as ~/libc/$ARCH/asa239.o"
