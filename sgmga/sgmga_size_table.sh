#!/bin/bash
#
gcc -c -g -I/$HOME/include sgmga_size_table.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_size_table.c."
  exit
fi
rm compiler.txt
#
gcc sgmga_size_table.o /$HOME/libc/$ARCH/sgmga.o /$HOME/libc/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_size_table.o."
  exit
fi
#
rm sgmga_size_table.o
#
mv a.out sgmga_size_table
./sgmga_size_table > sgmga_size_table_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_size_table."
  exit
fi
rm sgmga_size_table
#
echo "Program output written to sgmga_size_table_output.txt"
