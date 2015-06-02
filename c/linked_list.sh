#!/bin/bash
#
gcc -c linked_list.c
if [ $? -ne 0 ]; then
  echo "Errors compiling linked_list.c."
  exit
fi
#
gcc linked_list.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking linked_list.o."
  exit
fi
#
rm linked_list.o
#
mv a.out linked_list
./linked_list > linked_list_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running linked_list."
  exit
fi
rm linked_list
#
echo "Program output written to linked_list_output.txt"
