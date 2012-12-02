#!/bin/bash
#
cp table_io.h /$HOME/include
#
gcc -c -g table_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_io.c."
  exit
fi
rm compiler.txt
#
mv table_io.o ~/libc/$ARCH/table_io.o
#
echo "Library installed as ~/libc/$ARCH/table_io.o"
