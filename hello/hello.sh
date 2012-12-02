#!/bin/bash
#
gcc -c hello.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hello.c."
  exit
fi
rm compiler.txt
#
gcc hello.o
if [ $? -ne 0 ]; then
  echo "Errors linking hello.o."
  exit
fi
#
rm hello.o
#
mv a.out hello
./hello > hello_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hello."
  exit
fi
rm hello
#
echo "Program output written to hello_output.txt"
