#!/bin/bash
#
gcc -c -g mutex1.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mutex1.c"
  exit
fi
rm compiler.txt
#
gcc mutex1.o
if [ $? -ne 0 ]; then
  echo "Errors loading mutex1.c"
  exit
fi
#
rm mutex1.o
mv a.out mutex1
mutex1 > mutex1_output.txt
rm mutex1
#
echo "Program output written to mutex1_output.txt"
