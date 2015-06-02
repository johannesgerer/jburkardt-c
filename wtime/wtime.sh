#!/bin/bash
#
cp wtime.h /$HOME/include
#
gcc -c wtime.c
if [ $? -ne 0 ]; then
  echo "Errors compiling wtime.c."
  exit
fi
#
mv wtime.o ~/libc/$ARCH/wtime.o
#
echo "Library installed as ~/libc/$ARCH/wtime.o"
