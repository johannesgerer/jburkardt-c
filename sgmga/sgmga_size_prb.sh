#!/bin/bash
#
gcc -c -g -I/$HOME/include sgmga_size_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_size_prb.c."
  exit
fi
rm compiler.txt
#
gcc sgmga_size_prb.o /$HOME/libc/$ARCH/sgmga.o /$HOME/libc/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_size_prb.o."
  exit
fi
#
rm sgmga_size_prb.o
#
mv a.out sgmga_size_prb
./sgmga_size_prb > sgmga_size_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_size_prb."
  exit
fi
rm sgmga_size_prb
#
echo "Program output written to sgmga_size_prb_output.txt"
