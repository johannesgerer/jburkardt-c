#!/bin/bash
#
gcc -c -g ice_io_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ice_io_prb.c."
  exit
fi
rm compiler.txt
#
gcc ice_io_prb.o /$HOME/libc/$ARCH/ice_io.o -lnetcdf -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ice_io_prb.o."
  exit
fi
#
rm ice_io_prb.o
#
mv a.out ice_io_prb
./ice_io_prb > ice_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ice_io_prb."
  exit
fi
rm ice_io_prb
#
echo "Program output written to ice_io_prb_output.txt"
