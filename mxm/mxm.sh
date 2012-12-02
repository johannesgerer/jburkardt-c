#!/bin/bash
#
gcc -c -g mxm.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mxm.c."
  exit
fi
rm compiler.txt
#
gcc mxm.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mxm.o."
  exit
fi
#
rm mxm.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/mxm
#
echo "Executable installed as ~/binc/$ARCH/mxm"
