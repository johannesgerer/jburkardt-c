#!/bin/bash
#
cp hypercube_integrals.h /$HOME/include
#
gcc -c -g -I /$HOME/include hypercube_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hypercube_integrals.c"
  exit
fi
rm compiler.txt
#
mv hypercube_integrals.o ~/libc/$ARCH/hypercube_integrals.o
#
echo "Library installed as ~/libc/$ARCH/hypercube_integrals.o"
