#!/bin/bash
#
gcc -c sizes.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sizes.c"
  exit
fi
rm compiler.txt
#
gcc sizes.o
if [ $? -ne 0 ]; then
  echo "Errors loading sizes.o"
  exit
fi
rm sizes.o
#
mv a.out sizes
./sizes > sizes_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sizes"
  exit
fi
rm sizes
#
echo "Program output written to sizes_output.txt"
