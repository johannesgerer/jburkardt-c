#!/bin/bash
#
gcc -c -g -I/$HOME/include mesh_io_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mesh_io_prb.c."
  exit
fi
rm compiler.txt
#
gcc mesh_io_prb.o /$HOME/libc/$ARCH/mesh_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mesh_io_prb.o."
  exit
fi
#
rm mesh_io_prb.o
#
mv a.out mesh_io_prb
./mesh_io_prb > mesh_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running mesh_io_prb."
  exit
fi
rm mesh_io_prb
#
echo "Program output written to mesh_io_prb_output.txt"
