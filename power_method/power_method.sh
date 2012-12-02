#!/bin/bash
#
cp power_method.h /$HOME/include
#
gcc -c -g power_method.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling power_method.c."
  exit
fi
rm compiler.txt
#
mv power_method.o ~/libc/$ARCH/power_method.o
#
echo "Library installed as ~/libc/$ARCH/power_method.o"
