#!/bin/bash
#
gcc -c linked_list_vector.c
if [ $? -ne 0 ]; then
  echo "Errors compiling linked_list_vector.c."
  exit
fi
#
gcc linked_list_vector.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking linked_list_vector.o."
  exit
fi
#
rm linked_list_vector.o
#
mv a.out linked_list_vector
./linked_list_vector > linked_list_vector_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running linked_list_vector."
  exit
fi
rm linked_list_vector
#
echo "Program output written to linked_list_vector_output.txt"
