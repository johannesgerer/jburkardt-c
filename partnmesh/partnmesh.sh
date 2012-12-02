#!/bin/bash
#
#  Copy the INCLUDE files.
#
cp ../metis/*.h .
#
gcc -c -g -I . partnmesh.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling partnmesh.c."
  exit
fi
rm compiler.txt
#
gcc partnmesh.o -L$HOME/libc/$ARCH -lmetis -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading the object files."
  exit
fi
#
rm *.h
rm *.o
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/partnmesh
#
echo "Executable installed as ~/binc/$ARCH/partnmesh"
