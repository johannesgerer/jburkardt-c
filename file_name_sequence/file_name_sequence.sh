#!/bin/bash
#
cp file_name_sequence.h /$HOME/include
#
gcc -c -g file_name_sequence.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling file_name_sequence.c."
  exit
fi
rm compiler.txt
#
mv file_name_sequence.o ~/libc/$ARCH/file_name_sequence.o
#
echo "Library installed as ~/libc/$ARCH/file_name_sequence.o"
