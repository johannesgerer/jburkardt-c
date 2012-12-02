#!/bin/bash
#
gcc -c f77split.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling f77split.c."
  exit
fi
rm compiler.txt
#
gcc f77split.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading f77split.o."
  exit
fi
#
rm f77split.o
#
chmod u+x a.out
mv a.out ~/binc/$ARCH/f77split
#
echo "Executable installed as ~/binc/$ARCH/f77split"
