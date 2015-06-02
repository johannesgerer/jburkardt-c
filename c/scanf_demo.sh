#!/bin/bash
#
gcc -c scanf_demo.c
if [ $? -ne 0 ]; then
  echo "Errors compiling scanf_demo.c"
  exit
fi
#
gcc scanf_demo.o
if [ $? -ne 0 ]; then
  echo "Errors loading scanf_demo.o"
  exit
fi
rm scanf_demo.o
#
mv a.out scanf_demo
./scanf_demo < scanf_demo_input.txt > scanf_demo_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running scanf_demo"
  exit
fi
rm scanf_demo
#
echo "Program output written to scanf_demo_output.txt"
