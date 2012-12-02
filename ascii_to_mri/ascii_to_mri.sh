#!/bin/bash
#
gcc -c -g ascii_to_mri.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ascii_to_mri.c"
  exit
fi
rm compiler.txt
#
gcc ascii_to_mri.o
if [ $? -ne 0 ]; then
  echo "Errors loading ascii_to_mri.c"
  exit
fi
#
rm ascii_to_mri.o
mv a.out ~/binc/$ARCH/ascii_to_mri
#
echo "Executable installed as ~/binc/$ARCH/ascii_to_mri"
