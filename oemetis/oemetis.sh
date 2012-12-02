#!/bin/bash
#
#  Copy the INCLUDE files.
#
cp ../metis/*.h .
#
gcc -c -g -I . oemetis.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling oemetis.c."
  exit
fi
rm compiler.txt
#
gcc oemetis.o -L$HOME/libc/$ARCH -lmetis -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading the object files."
  exit
fi
#
rm *.h
rm *.o
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/oemetis
#
echo "Executable installed as ~/binc/$ARCH/oemetis"
