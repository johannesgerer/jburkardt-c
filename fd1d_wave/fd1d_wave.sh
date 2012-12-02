#!/bin/bash
#
cp fd1d_wave.h /$HOME/include
#
gcc -c -g -I /$HOME/include fd1d_wave.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_wave.c"
  exit
fi
rm compiler.txt
#
mv fd1d_wave.o ~/libc/$ARCH/fd1d_wave.o
#
echo "Library installed as ~/libc/$ARCH/fd1d_wave.o"
