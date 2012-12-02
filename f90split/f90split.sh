#!/bin/bash
#
gcc -c f90split.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling f90split.c."
  exit
fi
rm compiler.txt
#
gcc f90split.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading f90split.o."
  exit
fi
#
rm f90split.o
#
chmod u+x a.out
mv a.out ~/binc/$ARCH/f90split
#
echo "Executable installed as ~/binc/$ARCH/f90split"
