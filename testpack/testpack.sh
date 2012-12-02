#!/bin/bash
#
gcc -c -g testpack.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling testpack.c"
  exit
fi
rm compiler.txt
#
gcc testpack.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading testpack.o."
  exit
fi
#
rm testpack.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/testpack
#
echo "Executable installed as ~/binc/$ARCH/testpack"
