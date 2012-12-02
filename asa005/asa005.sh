#!/bin/bash
#
cp asa005.h /$HOME/include
#
gcc -c -g asa005.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa005.c."
  exit
fi
rm compiler.txt
#
mv asa005.o ~/libc/$ARCH/asa005.o
#
echo "Library installed as ~/libc/$ARCH/asa005.o"
