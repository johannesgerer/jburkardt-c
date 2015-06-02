#!/bin/bash
#
cp svd_snowfall.h /$HOME/include
#
gcc -c -g -I/$HOME/include svd_snowfall.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling svd_snowfall.c"
  exit
fi
rm compiler.txt
#
mv svd_snowfall.o ~/libc/$ARCH/svd_snowfall.o
#
echo "Library installed as ~/libc/$ARCH/svd_snowfall.o"
