#!/bin/bash
#
cp square_exactness.h /$HOME/include
#
gcc -c -I/$HOME/include square_exactness.c
if [ $? -ne 0 ]; then
  echo "Errors compiling square_exactness.c"
  exit
fi
#
mv square_exactness.o ~/libc/$ARCH/square_exactness.o
#
echo "Library installed as ~/libc/$ARCH/square_exactness.o"
