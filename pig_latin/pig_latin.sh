#!/bin/bash
#
gcc -c -g pig_latin.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pig_latin.c"
  exit
fi
rm compiler.txt
#
gcc pig_latin.o
if [ $? -ne 0 ]; then
  echo "Errors loading pig_latin.c"
  exit
fi
#
rm pig_latin.o
mv a.out ~/binc/$ARCH/pig_latin
#
echo "Executable installed as ~/binc/$ARCH/pig_latin"
