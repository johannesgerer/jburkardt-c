#!/bin/bash
#
gcc -c character_arithmetic.c
if [ $? -ne 0 ]; then
  echo "Errors compiling character_arithmetic.c."
  exit
fi
#
gcc character_arithmetic.o
if [ $? -ne 0 ]; then
  echo "Errors linking character_arithmetic.o."
  exit
fi
#
rm character_arithmetic.o
#
mv a.out character_arithmetic
./character_arithmetic > character_arithmetic_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running character_arithmetic."
  exit
fi
rm character_arithmetic
#
echo "Program output written to character_arithmetic_output.txt"
