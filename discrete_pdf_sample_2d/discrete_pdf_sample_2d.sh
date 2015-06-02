#!/bin/bash
#
gcc -c -g discrete_pdf_sample_2d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling discrete_pdf_sample_2d.c"
  exit
fi
rm compiler.txt
#
gcc discrete_pdf_sample_2d.o
if [ $? -ne 0 ]; then
  echo "Errors while loading discrete_pdf_sample_2d.o"
  exit
fi
rm discrete_pdf_sample_2d.o
#
mv a.out ~/binc/$ARCH/discrete_pdf_sample_2d
#
echo "Executable installed as ~/binc/$ARCH/discrete_pdf_sample_2d"
