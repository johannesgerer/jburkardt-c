#!/bin/bash
#
cp table_io.h /$HOME/include
#
gcc -c table_io.c
if [ $? -ne 0 ]; then
  echo "Errors compiling table_io.c."
  exit
fi
#
mv table_io.o ~/libc/$ARCH/table_io.o
#
echo "Library installed as ~/libc/$ARCH/table_io.o"
