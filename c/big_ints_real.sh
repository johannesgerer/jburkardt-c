#!/bin/bash
#
gcc -c big_ints_real.c
if [ $? -ne 0 ]; then
  echo "Errors compiling big_ints_real.c."
  exit
fi
#
gcc big_ints_real.o
if [ $? -ne 0 ]; then
  echo "Errors linking big_ints_real.o."
  exit
fi
#
rm big_ints_real.o
#
mv a.out big_ints_real
./big_ints_real > big_ints_real_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running big_ints_real."
  exit
fi
rm big_ints_real
#
echo "Program output written to big_ints_real_output.txt"
