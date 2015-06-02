#!/bin/bash
#
gcc -c ziggurat_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling ziggurat_prb.c."
  exit
fi
#
gcc ziggurat_prb.o /$HOME/libc/$ARCH/ziggurat.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ziggurat_prb.o."
  exit
fi
#
rm ziggurat_prb.o
#
mv a.out ziggurat_prb
./ziggurat_prb > ziggurat_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ziggurat_prb."
  exit
fi
rm ziggurat_prb
#
echo "Program output written to ziggurat_prb_output.txt"
