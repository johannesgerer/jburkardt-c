#!/bin/bash
#
cp bvec.h /$HOME/include
#
gcc -c bvec.c
if [ $? -ne 0 ]; then
  echo "Errors compiling bvec.c."
  exit
fi
#
mv bvec.o ~/libc/$ARCH/bvec.o
#
echo "Library installed as ~/libc/$ARCH/bvec.o"
