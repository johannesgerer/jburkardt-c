#!/bin/bash
#
cp monomial_value.h /$HOME/include
#
gcc -c -I/$HOME/include monomial_value.c
if [ $? -ne 0 ]; then
  echo "Errors compiling monomial_value.c"
  exit
fi
#
mv monomial_value.o ~/libc/$ARCH/monomial_value.o
#
echo "Library installed as ~/libc/$ARCH/monomial_value.o"
