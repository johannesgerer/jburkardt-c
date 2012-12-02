#!/bin/bash
#
gcc -c -g nas.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nas.c"
  exit
fi
rm compiler.txt
#
gcc nas.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nas.o"
  exit
fi
rm nas.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/nas
#
echo "Executable installed as ~/binc/$ARCH/nas"
