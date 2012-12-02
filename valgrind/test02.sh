#!/bin/bash
#
gcc -c -g test02.c
if [ $? -ne 0 ]; then
  echo "Errors compiling test02.c"
  exit
fi
#
gcc test02.o
if [ $? -ne 0 ]; then
  echo "Errors loading test02.c"
  exit
fi
#
mv a.out test02
echo "Executable created as test02"
#
#  Run program
#
test02 > test02_output.txt
echo "Executable output stored as test02_output.txt"
#
#  Rerun program with valgrind
#
valgrind --dsymutil=yes test02 &> test02_valgrind_output.txt
echo "Valgrind output stored as test02_valgrind_output.txt"
#
#  Don't delete object code til later!
#
rm test02.o
