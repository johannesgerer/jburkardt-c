#!/bin/bash
#
cp disk_integrals.h /$HOME/include
#
gcc -c -g -I/$HOME/include disk_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling disk_integrals.c"
  exit
fi
rm compiler.txt
#
mv disk_integrals.o ~/libc/$ARCH/disk_integrals.o
#
echo "Library installed as ~/libc/$ARCH/disk_integrals.o"
