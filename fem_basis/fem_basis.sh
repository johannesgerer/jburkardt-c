#!/bin/bash
#
cp fem_basis.h /$HOME/include
#
gcc -c -I /$HOME/include fem_basis.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_basis.c"
  exit
fi
#
mv fem_basis.o ~/libc/$ARCH/fem_basis.o
#
echo "Library installed as ~/libc/$ARCH/fem_basis.o"
