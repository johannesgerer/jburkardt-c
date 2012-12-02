#!/bin/bash
#
cp test_ls.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_ls.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_ls.c"
  exit
fi
rm compiler.txt
#
mv test_ls.o ~/libc/$ARCH/test_ls.o
#
echo "Library installed as ~/libc/$ARCH/test_ls.o"
