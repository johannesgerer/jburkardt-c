#!/bin/bash
#
gcc -c linpack_bench.c
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_bench.c."
  exit
fi
#
gcc linpack_bench.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linpack_bench.o."
  exit
fi
#
rm linpack_bench.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/linpack_bench
#
echo "Executable installed as ~/binc/$ARCH/linpack_bench"
