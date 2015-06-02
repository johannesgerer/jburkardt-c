#!/bin/bash
#
gcc -c function_pointer.c
if [ $? -ne 0 ]; then
  echo "Errors while compiling function_pointer.c"
  exit
fi
#
gcc function_pointer.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading function_pointer.o"
  exit
fi
rm function_pointer.o
#
mv a.out function_pointer
./function_pointer > function_pointer_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running function_pointer"
  exit
fi
rm function_pointer
#
echo "Program output written to function_pointer_output.txt"
