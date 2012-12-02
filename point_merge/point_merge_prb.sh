#!/bin/bash
#
gcc -c -g point_merge_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling point_merge_prb.c."
  exit
fi
rm compiler.txt
#
gcc point_merge_prb.o /$HOME/libc/$ARCH/point_merge.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading point_merge_prb.o."
  exit
fi
#
rm point_merge_prb.o
#
mv a.out point_merge_prb
./point_merge_prb > point_merge_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running point_merge_prb."
  exit
fi
rm point_merge_prb
#
echo "Program output written to point_merge_prb_output.txt"
