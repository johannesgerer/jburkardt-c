#!/bin/bash
#
cp asa266.h /$HOME/include
#
gcc -c -g asa266.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa266.c."
  exit
fi
rm compiler.txt
#
mv asa266.o ~/libc/$ARCH/asa266.o
#
echo "Library installed as ~/libc/$ARCH/asa266.o"
