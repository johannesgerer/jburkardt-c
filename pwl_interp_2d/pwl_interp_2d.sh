#!/bin/bash
#
cp pwl_interp_2d.h /$HOME/include
#
gcc -c -g -I/$HOME/include pwl_interp_2d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_interp_2d.c"
  exit
fi
rm compiler.txt
#
mv pwl_interp_2d.o ~/libc/$ARCH/pwl_interp_2d.o
#
echo "Library installed as ~/libc/$ARCH/pwl_interp_2d.o"
