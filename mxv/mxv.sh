#!/bin/bash
#
gcc -c mxv.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mxv.c"
  exit
fi
rm compiler.txt
#
gcc mxv.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mxv.o"
  exit
fi
rm mxv.o
#
mv a.out ~/binc/$ARCH/mxv
#
echo "Executable installed as ~/binc/$ARCH/mxv"
