#!/bin/bash
#
gcc -c -g -I/$HOME/include problem1.c
if [ $? -ne 0 ]; then
  echo "Errors compiling problem1.c"
  exit
fi
#
gcc -c -g -I/$HOME/include problem1_main.c
if [ $? -ne 0 ]; then
  echo "Errors compiling problem1_main.c"
  exit
fi
#
gcc problem1_main.o \
  problem1.o \
  /$HOME/libc/$ARCH/pdflib.o \
  /$HOME/libc/$ARCH/rnglib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading problem1_main.o + problem1.o + pdflib + rnglib"
  exit
fi
rm problem1_main.o
rm problem1.o
#
mv a.out problem1_main
./problem1_main > problem1_main_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running problem1_main"
  exit
fi
rm problem1_main
#
echo "Test program output written to problem1_main_output.txt."
