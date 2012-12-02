#!/bin/bash
#
cp black_scholes.h /$HOME/include
#
gcc -c -g -I /$HOME/include black_scholes.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling black_scholes.c"
  exit
fi
rm compiler.txt
#
mv black_scholes.o ~/libc/$ARCH/black_scholes.o
#
echo "Library installed as ~/libc/$ARCH/black_scholes.o"
