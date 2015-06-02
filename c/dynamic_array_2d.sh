#!/bin/bash
#
gcc -c dynamic_array_2d.c
if [ $? -ne 0 ]; then
  echo "Errors compiling dynamic_array_2d.c."
  exit
fi
#
gcc dynamic_array_2d.o
if [ $? -ne 0 ]; then
  echo "Errors linking dynamic_array_2d.o."
  exit
fi
#
rm dynamic_array_2d.o
#
mv a.out dynamic_array_2d
./dynamic_array_2d > dynamic_array_2d_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dynamic_array_2d."
  exit
fi
rm dynamic_array_2d
#
echo "Program output written to dynamic_array_2d_output.txt"
