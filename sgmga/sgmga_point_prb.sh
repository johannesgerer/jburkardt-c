#!/bin/bash
#
gcc -c -g -I/$HOME/include sgmga_point_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_point_prb.c."
  exit
fi
rm compiler.txt
#
gcc sgmga_point_prb.o /$HOME/libc/$ARCH/sgmga.o /$HOME/libc/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_point_prb.o."
  exit
fi
#
rm sgmga_point_prb.o
#
mv a.out sgmga_point_prb
./sgmga_point_prb > sgmga_point_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_point_prb."
  exit
fi
rm sgmga_point_prb
#
echo "Program output written to sgmga_point_prb_output.txt"
