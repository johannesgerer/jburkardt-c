#!/bin/bash
#
cp fem1d_lagrange.h /$HOME/include
#
gcc -c -I/$HOME/include fem1d_lagrange.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_lagrange.c"
  exit
fi
#
mv fem1d_lagrange.o ~/libc/$ARCH/fem1d_lagrange.o
#
echo "Library installed as ~/libc/$ARCH/fem1d_lagrange.o"
