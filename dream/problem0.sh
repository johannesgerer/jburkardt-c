#!/bin/bash
#
gcc -c -g -I/$HOME/include problem0.c
if [ $? -ne 0 ]; then
  echo "Errors compiling problem0.c"
  exit
fi
#
gcc /$HOME/libc/$ARCH/dream.o \
  problem0.o \
  /$HOME/libc/$ARCH/pdflib.o \
  /$HOME/libc/$ARCH/rnglib.o -lm
#
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dream.o + problem0.o + pdflib.o + rnglib.o"
  exit
fi
#
rm problem0.o
#
mv a.out problem0
./problem0 > problem0_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running problem0."
  exit
fi
rm problem0
#
echo "Program output written to problem0_output.txt"
