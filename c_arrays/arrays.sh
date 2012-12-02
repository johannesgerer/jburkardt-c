#!/bin/bash
#
gcc -c arrays.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling arrays.c."
  exit
fi
rm compiler.txt
#
gcc arrays.o
if [ $? -ne 0 ]; then
  echo "Errors linking arrays.o."
  exit
fi
#
rm arrays.o
#
mv a.out arrays
./arrays > arrays_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running arrays."
  exit
fi
rm arrays
#
echo "Program output written to arrays_output.txt"
