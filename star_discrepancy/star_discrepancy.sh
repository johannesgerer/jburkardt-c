#!/bin/bash
#
gcc -c star_discrepancy.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling star_discrepancy.c"
  exit
fi
rm compiler.txt
#
gcc star_discrepancy.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading star_discrepancy.o."
  exit
fi
#
rm star_discrepancy.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/star_discrepancy
#
echo "Executable installed as ~/binc/$ARCH/star_discrepancy"
