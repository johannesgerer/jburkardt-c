#!/bin/bash
#
cp hypersphere_properties.h /$HOME/include
#
gcc -c -g -I/$HOME/include hypersphere_properties.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hypersphere_properties.c"
  exit
fi
rm compiler.txt
#
mv hypersphere_properties.o ~/libc/$ARCH/hypersphere_properties.o
#
echo "Library installed as ~/libc/$ARCH/hypersphere_properties.o"
