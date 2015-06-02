#!/bin/bash
#
gcc -c fem_to_medit.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_to_medit.c"
  exit
fi
#
gcc fem_to_medit.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_to_medit.o"
  exit
fi
#
rm fem_to_medit.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/fem_to_medit
#
echo "Program installed as ~/binc/$ARCH/fem_to_medit"
