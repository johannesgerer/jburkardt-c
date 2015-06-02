#!/bin/bash
#
gcc -c -I/$HOME/include problem1.c
if [ $? -ne 0 ]; then
  echo "Errors compiling problem1.c"
  exit
fi
#
gcc /$HOME/libc/$ARCH/dream.o \
  problem1.o \
  /$HOME/libc/$ARCH/pdflib.o \
  /$HOME/libc/$ARCH/rnglib.o -lm
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dream.o + problem1.o + pdflib.o + rnglib.o"
  exit
fi
#
rm problem1.o
#
mv a.out problem1
./problem1 > problem1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running problem1."
  exit
fi
rm problem1
#
echo "Program output written to problem1_output.txt"
