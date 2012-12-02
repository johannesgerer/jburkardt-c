#!/bin/bash
#
cp asa066.h /$HOME/include
#
gcc -c -g asa066.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa066.c."
  exit
fi
rm compiler.txt
#
mv asa066.o ~/libc/$ARCH/asa066.o
#
echo "Library installed as ~/libc/$ARCH/asa066.o"
