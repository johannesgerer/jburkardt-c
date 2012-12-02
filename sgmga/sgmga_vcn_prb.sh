#!/bin/bash
#
gcc -c -g -I/$HOME/include sgmga_vcn_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_vcn_prb.c."
  exit
fi
rm compiler.txt
#
gcc sgmga_vcn_prb.o /$HOME/libc/$ARCH/sgmga.o /$HOME/libc/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_vcn_prb.o."
  exit
fi
#
rm sgmga_vcn_prb.o
#
mv a.out sgmga_vcn_prb
./sgmga_vcn_prb > sgmga_vcn_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_vcn_prb."
  exit
fi
rm sgmga_vcn_prb
#
echo "Program output written to sgmga_vcn_prb_output.txt"
