#!/bin/bash
#
gcc -c -g -I/$HOME/include ppmb_io_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ppmb_io_prb.c."
  exit
fi
rm compiler.txt
#
gcc ppmb_io_prb.o /$HOME/libc/$ARCH/ppmb_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ppmb_io_prb.o."
  exit
fi
#
rm ppmb_io_prb.o
#
mv a.out ppmb_io_prb
./ppmb_io_prb > ppmb_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ppmb_io_prb."
  exit
fi
rm ppmb_io_prb
#
echo "Program output written to ppmb_io_prb_output.txt"
