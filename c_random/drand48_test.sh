#!/bin/bash
#
gcc -c drand48_test.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling drand48_test.c"
  exit
fi
rm compiler.txt
#
gcc drand48_test.o
if [ $? -ne 0 ]; then
  echo "Errors loading drand48_test.o"
  exit
fi
rm drand48_test.o
#
mv a.out drand48_test
./drand48_test > drand48_test_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running drand48_test"
  exit
fi
rm drand48_test
#
echo "Program output written to drand48_test_output.txt"
