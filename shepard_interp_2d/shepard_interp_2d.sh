#!/bin/bash
#
cp shepard_interp_2d.h /$HOME/include
#
gcc -c -g -I/$HOME/include shepard_interp_2d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling shepard_interp_2d.c"
  exit
fi
rm compiler.txt
#
mv shepard_interp_2d.o ~/libc/$ARCH/shepard_interp_2d.o
#
echo "Library installed as ~/libc/$ARCH/shepard_interp_2d.o"
