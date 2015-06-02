#!/bin/bash
#
cp asa109.h /$HOME/include
#
gcc -c asa109.c
if [ $? -ne 0 ]; then
  echo "Errors compiling asa109.c."
  exit
fi
#
mv asa109.o ~/libc/$ARCH/asa109.o
#
echo "Library installed as ~/libc/$ARCH/asa109.o"
