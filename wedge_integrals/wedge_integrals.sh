#!/bin/bash
#
cp wedge_integrals.h /$HOME/include
#
gcc -c -I/$HOME/include wedge_integrals.c
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_integrals.c"
  exit
fi
#
mv wedge_integrals.o ~/libc/$ARCH/wedge_integrals.o
#
echo "Library installed as ~/libc/$ARCH/wedge_integrals.o"
