#!/bin/bash
#
gcc -c -g -I/$HOME/include wishart_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wishart_prb.c"
  exit
fi
rm compiler.txt
#
gcc wishart_prb.o \
  /$HOME/libc/$ARCH/wishart.o \
  /$HOME/libc/$ARCH/pdflib.o \
  /$HOME/libc/$ARCH/rnglib.o -lm

if [ $? -ne 0 ]; then
  echo "Errors linking and loading wishart_prb.o."
  exit
fi
#
rm wishart_prb.o
#
mv a.out wishart_prb
./wishart_prb > wishart_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wishart_prb."
  exit
fi
rm wishart_prb
#
echo "Program output written to wishart_prb_output.txt"
