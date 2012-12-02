#!/bin/bash
#
gcc -c -O2 md.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling md.c"
  exit
fi
rm compiler.txt
#
gcc md.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading md.o"
  exit
fi
rm md.o
#
mv a.out md_O2
./md_O2 < md_input.txt > md_O2_output.txt
rm md_O2
#
echo "Output written to md_O2_output.txt"
