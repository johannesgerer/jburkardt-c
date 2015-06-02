#!/bin/bash
#
gcc -c -I/$HOME/include fem_basis_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_basis_prb.c"
  exit
fi
#
gcc fem_basis_prb.o /$HOME/libc/$ARCH/fem_basis.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_basis_prb.o."
  exit
fi
#
rm fem_basis_prb.o
#
mv a.out fem_basis_prb
./fem_basis_prb > fem_basis_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem_basis_prb."
  exit
fi
rm fem_basis_prb
#
echo "Program output written to fem_basis_prb_output.txt"
