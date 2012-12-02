#!/bin/bash
#
#  Copy the INCLUDE files.
#
cp ../metis/*.h .
#
gcc -c -g -I . mesh2nodal.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mesh2nodal.c."
  exit
fi
rm compiler.txt
#
gcc mesh2nodal.o -L$HOME/libc/$ARCH -lmetis -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading the object files."
  exit
fi
#
rm *.h
rm *.o
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/mesh2nodal
#
echo "Executable installed as ~/binc/$ARCH/mesh2nodal"
