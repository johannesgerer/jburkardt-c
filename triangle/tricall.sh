#!/bin/bash
#
gcc -c tricall.c
if [ $? -ne 0 ]; then
  echo "Errors compiling tricall.c."
  exit
fi
#
gcc tricall.o /$HOME/libc/$ARCH/triangle.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tricall.o."
  exit
fi
#
rm tricall.o
#
mv a.out tricall
./tricall > tricall_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tricall."
  exit
fi
rm tricall
#
echo "Program output written to tricall_output.txt"
