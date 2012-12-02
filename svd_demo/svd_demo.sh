#!/bin/bash
#
gcc -c -g -I$HOME/include svd_demo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling svd_demo.c"
  exit
fi
rm compiler.txt
#
gcc svd_demo.o $HOME/libc/$ARCH/linpack_d.o $HOME/libc/$ARCH/blas1_d.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading svd_demo.o + linpack_d.o + blas1_d.o."
  exit
fi
#
rm svd_demo.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/svd_demo
#
echo "Executable installed as ~/binc/$ARCH/svd_demo"
