#!/bin/bash
#
cp polygon_properties.h /$HOME/include
#
gcc -c -I/$HOME/include polygon_properties.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_properties.c"
  exit
fi
rm compiler.txt
#
mv polygon_properties.o ~/libc/$ARCH/polygon_properties.o
#
echo "Library installed as ~/libc/$ARCH/polygon_properties.o"
