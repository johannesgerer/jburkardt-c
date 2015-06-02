#!/bin/bash
#
gcc -c laguerre_exactness.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laguerre_exactness.c"
  exit
fi
rm compiler.txt
#
gcc laguerre_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading laguerre_exactness.o"
  exit
fi
rm laguerre_exactness.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/laguerre_exactness
#
echo "Executable installed as ~/binc/$ARCH/laguerre_exactness"
