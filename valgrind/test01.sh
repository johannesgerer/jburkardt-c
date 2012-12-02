#!/bin/bash
#
gcc -c -g test01.c
if [ $? -ne 0 ]; then
  echo "Errors compiling test01.c"
  exit
fi
#
gcc test01.o
if [ $? -ne 0 ]; then
  echo "Errors loading test01.c"
  exit
fi
#
mv a.out test01
echo "Executable created as test01"
#
#  Run program
#
test01 > test01_output.txt
echo "Executable output stored as test01_output.txt"
#
#  Rerun program with valgrind
#
valgrind --dsymutil=yes test01 &> test01_valgrind_output.txt
echo "Valgrind output stored as test01_valgrind_output.txt"
#
#  Don't delete object code til later!
#
rm test01.o
