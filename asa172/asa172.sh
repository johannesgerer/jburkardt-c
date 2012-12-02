#!/bin/bash
#
cp asa172.h /$HOME/include
#
gcc -c -g asa172.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa172.c."
  exit
fi
rm compiler.txt
#
mv asa172.o ~/libc/$ARCH/asa172.o
#
echo "Library installed as ~/libc/$ARCH/asa172.o"
