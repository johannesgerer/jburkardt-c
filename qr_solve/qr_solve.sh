#!/bin/bash
#
cp qr_solve.h /$HOME/include
#
gcc -c -g -I /$HOME/include qr_solve.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qr_solve.c"
  exit
fi
rm compiler.txt
#
mv qr_solve.o ~/libc/$ARCH/qr_solve.o
#
echo "Library installed as ~/libc/$ARCH/qr_solve.o"
