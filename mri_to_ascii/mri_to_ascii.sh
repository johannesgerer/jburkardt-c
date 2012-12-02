#!/bin/bash
#
gcc -c -g mri_to_ascii.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mri_to_ascii.c"
  exit
fi
rm compiler.txt
#
gcc mri_to_ascii.o
if [ $? -ne 0 ]; then
  echo "Errors loading mri_to_ascii.c"
  exit
fi
#
rm mri_to_ascii.o
mv a.out ~/binc/$ARCH/mri_to_ascii
#
echo "Executable installed as ~/binc/$ARCH/mri_to_ascii"
