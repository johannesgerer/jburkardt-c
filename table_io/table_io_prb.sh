#!/bin/bash
#
gcc -c table_io_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling table_io_prb.c."
  exit
fi
#
gcc table_io_prb.o /$HOME/libc/$ARCH/table_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_io_prb.o."
  exit
fi
#
rm table_io_prb.o
#
mv a.out table_io_prb
./table_io_prb > table_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running table_io_prb."
  exit
fi
rm table_io_prb
#
echo "Program output written to table_io_prb_output.txt"
