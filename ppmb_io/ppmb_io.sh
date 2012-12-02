#!/bin/bash
#
cp ppmb_io.h /$HOME/include
#
gcc -c -g -I /$HOME/include ppmb_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ppmb_io.c."
  exit
fi
rm compiler.txt
#
mv ppmb_io.o ~/libc/$ARCH/ppmb_io.o
#
echo "Library installed as ~/libc/$ARCH/ppmb_io.o"
