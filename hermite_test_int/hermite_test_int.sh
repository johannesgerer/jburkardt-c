#!/bin/bash
#
cp hermite_test_int.h /$HOME/include
#
gcc -c -I /$HOME/include hermite_test_int.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_test_int.c"
  exit
fi
rm compiler.txt
#
mv hermite_test_int.o ~/libc/$ARCH/hermite_test_int.o
#
echo "Library installed as ~/libc/$ARCH/hermite_test_int.o"
