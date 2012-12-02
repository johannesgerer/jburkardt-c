#!/bin/bash
#
#  Get INCLUDE files.
#
cp ../metis/*.h .
#
gcc -c -g -I . kmetis.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kmetis.c."
  exit
fi
rm compiler.txt
#
gcc kmetis.o -L$HOME/libc/$ARCH -lmetis -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading the object files."
  exit
fi
#
rm *.h
rm *.o
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/kmetis
#
echo "Executable installed as ~/binc/$ARCH/kmetis"
