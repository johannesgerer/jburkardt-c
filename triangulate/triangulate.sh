#!/bin/bash
#
gcc -c triangulate.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulate.c."
  exit
fi
rm compiler.txt
#
gcc triangulate.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulate.o."
  exit
fi
#
rm triangulate.o
#
chmod u+x a.out
mv a.out ~/binc/$ARCH/triangulate
#
echo "Executable installed as ~/binc/$ARCH/triangulate"
