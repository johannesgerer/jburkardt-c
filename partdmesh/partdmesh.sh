#!/bin/bash
#
#  Copy the INCLUDE files.
#
cp ../metis/*.h .
#
gcc -c -g -I . partdmesh.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling partdmesh.c."
  exit
fi
rm compiler.txt
#
gcc partdmesh.o -L$HOME/libc/$ARCH -lmetis -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading the object files."
  exit
fi
#
rm *.h
rm *.o
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/partdmesh
#
echo "Executable installed as ~/binc/$ARCH/partdmesh"
