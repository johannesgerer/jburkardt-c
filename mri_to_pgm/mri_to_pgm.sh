#!/bin/bash
#
gcc -c -g mri_to_pgm.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mri_to_pgm.c"
  exit
fi
rm compiler.txt
#
gcc mri_to_pgm.o
if [ $? -ne 0 ]; then
  echo "Errors loading mri_to_pgm.c"
  exit
fi
#
rm mri_to_pgm.o
mv a.out ~/binc/$ARCH/mri_to_pgm
#
echo "Executable installed as ~/binc/$ARCH/mri_to_pgm"
