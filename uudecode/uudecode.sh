#!/bin/bash
#
gcc -c -g uudecode.c
if [ $? -ne 0 ]; then
  echo "Errors compiling uudecode.c"
  exit
fi
#
gcc uudecode.o
if [ $? -ne 0 ]; then
  echo "Errors loading uudecode.c"
  exit
fi
#
rm uudecode.o
mv a.out ~/binc/$ARCH/uudecode
#
echo "Executable installed as ~/binc/$ARCH/uudecode."
