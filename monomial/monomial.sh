#!/bin/bash
#
cp monomial.h /$HOME/include
#
gcc -c -I/$HOME/include monomial.c
if [ $? -ne 0 ]; then
  echo "Errors compiling monomial.c"
  exit
fi
#
mv monomial.o ~/libc/$ARCH/monomial.o
#
echo "Library installed as ~/libc/$ARCH/monomial.o"
