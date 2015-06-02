#!/bin/bash
#
cp pdflib.h /$HOME/include
#
gcc -c -g -I/$HOME/include pdflib.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pdflib.c"
  exit
fi
rm compiler.txt
#
mv pdflib.o ~/libc/$ARCH/pdflib.o
#
echo "Library installed as ~/libc/$ARCH/pdflib.o"
