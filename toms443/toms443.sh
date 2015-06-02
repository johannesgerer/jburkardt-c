#!/bin/bash
#
cp toms443.h /$HOME/include
#
gcc -c -I/$HOME/include toms443.c
if [ $? -ne 0 ]; then
  echo "Errors compiling toms443.c"
  exit
fi
#
mv toms443.o ~/libc/$ARCH/toms443.o
#
echo "Library installed as ~/libc/$ARCH/toms443.o"
