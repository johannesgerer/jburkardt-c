#!/bin/bash
#
cp line_cvt_lloyd.h /$HOME/include
#
gcc -c -I/$HOME/include line_cvt_lloyd.c
if [ $? -ne 0 ]; then
  echo "Errors compiling line_cvt_lloyd.c"
  exit
fi
#
mv line_cvt_lloyd.o ~/libc/$ARCH/line_cvt_lloyd.o
#
echo "Library installed as ~/libc/$ARCH/line_cvt_lloyd.o"
