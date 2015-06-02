#!/bin/bash
#
cp latin_random.h /$HOME/include
#
gcc -c -I /$HOME/include latin_random.c
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_random.c"
  exit
fi
#
mv latin_random.o ~/libc/$ARCH/latin_random.o
#
echo "Library installed as ~/libc/$ARCH/latin_random.o"
