#!/bin/bash
#
cp prime_serial.h /$HOME/include
#
gcc -c -g prime_serial.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling prime_serial.c."
  exit
fi
rm compiler.txt
#
mv prime_serial.o ~/libc/$ARCH/prime_serial.o
#
echo "Library installed as ~/libc/$ARCH/prime_serial.o"
