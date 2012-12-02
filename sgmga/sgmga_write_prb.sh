#!/bin/bash
#
gcc -c -g -I/$HOME/include sgmga_write_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga_write_prb.c."
  exit
fi
rm compiler.txt
#
gcc sgmga_write_prb.o /$HOME/libc/$ARCH/sgmga.o /$HOME/libc/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sgmga_write_prb.o."
  exit
fi
#
rm sgmga_write_prb.o
#
mv a.out sgmga_write_prb
./sgmga_write_prb > sgmga_write_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sgmga_write_prb."
  exit
fi
rm sgmga_write_prb
#
echo "Program output written to sgmga_write_prb_output.txt"
#
#  Move sparse grid files to dataset directory.
#
mv *_a.txt ../../datasets/sgmga
mv *_n.txt ../../datasets/sgmga
mv *_p.txt ../../datasets/sgmga
mv *_r.txt ../../datasets/sgmga
mv *_w.txt ../../datasets/sgmga
mv *_x.txt ../../datasets/sgmga
#
echo "Program output files moved to ../../datasets/sgmga"
