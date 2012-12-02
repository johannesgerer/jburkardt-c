#!/bin/bash
#
cp brownian_motion_simulation.h /$HOME/include
#
gcc -c -g -I/$HOME/include brownian_motion_simulation.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling brownian_motion_simulation.c"
  exit
fi
rm compiler.txt
#
mv brownian_motion_simulation.o ~/libc/$ARCH/brownian_motion_simulation.o
#
echo "Library installed as ~/libc/$ARCH/brownian_motion_simulation.o"
