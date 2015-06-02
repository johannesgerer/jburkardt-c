#!/bin/bash
#
gcc -c -g -I/$HOME/include pdflib_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pdflib_prb.c"
  exit
fi
rm compiler.txt
#
gcc pdflib_prb.o /$HOME/libc/$ARCH/pdflib.o  /$HOME/libc/$ARCH/rnglib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pdflib_prb.o."
  exit
fi
#
rm pdflib_prb.o
#
mv a.out pdflib_prb
./pdflib_prb > pdflib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pdflib_prb."
  exit
fi
rm pdflib_prb
#
echo "Program output written to pdflib_prb_output.txt"
