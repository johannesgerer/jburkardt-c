#!/bin/bash
#
gcc -c -g -I/$HOME/include sgmga_index_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_index_prb.c."
  exit
fi
rm compiler.txt
#
gcc sgmga_index_prb.o /$HOME/libc/$ARCH/sgmga.o /$HOME/libc/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_index_prb.o."
  exit
fi
#
rm sgmga_index_prb.o
#
mv a.out sgmga_index_prb
./sgmga_index_prb > sgmga_index_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_index_prb."
  exit
fi
rm sgmga_index_prb
#
echo "Program output written to sgmga_index_prb_output.txt"
