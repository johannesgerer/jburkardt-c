#!/bin/bash
#
gcc -c -g uuencode.c
if [ $? -ne 0 ]; then
  echo "Errors compiling uuencode.c"
  exit
fi
#
gcc uuencode.o
if [ $? -ne 0 ]; then
  echo "Errors loading uuencode.c"
  exit
fi
#
rm uuencode.o
mv a.out ~/binc/$ARCH/uuencode
#
echo "Executable installed as ~/binc/$ARCH/uuencode."
