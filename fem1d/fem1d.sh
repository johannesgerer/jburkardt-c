#!/bin/bash
#
gcc -c -g fem1d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d.c."
  exit
fi
rm compiler.txt
#
gcc fem1d.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d.o."
  exit
fi
#
rm fem1d.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/fem1d
#
echo "Executable installed as ~/binc/$ARCH/fem1d"
