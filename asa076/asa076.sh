#!/bin/bash
#
cp asa076.h /$HOME/include
#
gcc -c -g asa076.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa076.c."
  exit
fi
rm compiler.txt
#
mv asa076.o ~/libc/$ARCH/asa076.o
#
echo "Library installed as ~/libc/$ARCH/asa076.o"
