#!/bin/bash
#
cp pwl_approx_1d.h /$HOME/include
#
gcc -c -g -I/$HOME/include pwl_approx_1d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_approx_1d.c"
  exit
fi
rm compiler.txt
#
mv pwl_approx_1d.o ~/libc/$ARCH/pwl_approx_1d.o
#
echo "Library installed as ~/libc/$ARCH/pwl_approx_1d.o"
