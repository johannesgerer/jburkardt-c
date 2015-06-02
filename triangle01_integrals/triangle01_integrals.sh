#!/bin/bash
#
cp triangle01_integrals.h /$HOME/include
#
gcc -c -I /$HOME/include triangle01_integrals.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle01_integrals.c"
  exit
fi
#
mv triangle01_integrals.o ~/libc/$ARCH/triangle01_integrals.o
#
echo "Library installed as ~/libc/$ARCH/triangle01_integrals.o"
