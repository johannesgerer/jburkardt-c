#!/bin/bash
#
cp jacobi_eigenvalue.h /$HOME/include
#
gcc -c -I/$HOME/include jacobi_eigenvalue.c
if [ $? -ne 0 ]; then
  echo "Errors compiling jacobi_eigenvalue.c"
  exit
fi
#
mv jacobi_eigenvalue.o ~/libc/$ARCH/jacobi_eigenvalue.o
#
echo "Library installed as ~/libc/$ARCH/jacobi_eigenvalue.o"
