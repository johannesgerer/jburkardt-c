#!/bin/bash
#
gcc -c -g feynman_kac_2d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling feynman_kac_2d.c"
  exit
fi
rm compiler.txt
#
gcc feynman_kac_2d.o
if [ $? -ne 0 ]; then
  echo "Errors loading feynman_kac_2d.c"
  exit
fi
#
rm feynman_kac_2d.o
mv a.out feynman_kac_2d
./feynman_kac_2d > feynman_kac_2d_output.txt
rm feynman_kac_2d
#
echo "Program output written to feynman_kac_2d_output.txt"
