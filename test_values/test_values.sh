#!/bin/bash
#
cp test_values.h /$HOME/include
#
gcc -c test_values.c
if [ $? -ne 0 ]; then
  echo "Errors compiling test_values.c."
  exit
fi
#
mv test_values.o ~/libc/$ARCH/test_values.o
#
echo "Library installed as ~/libc/$ARCH/test_values.o"
