#!/bin/bash
#
gcc -c sparse_grid_cc_dataset.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_cc_dataset.c."
  exit
fi
rm compiler.txt
#
gcc sparse_grid_cc_dataset.o
if [ $? -ne 0 ]; then
  echo "Errors loading sparse_grid_cc_dataset.c."
  exit
fi
mv a.out ~/binc/$ARCH/sparse_grid_cc_dataset
#
echo "Executable installed as ~/binc/$ARCH/sparse_grid_cc_dataset"
