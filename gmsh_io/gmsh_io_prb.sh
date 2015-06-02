#!/bin/bash
#
gcc -c -I/$HOME/include gmsh_io_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling gmsh_io_prb.c"
  exit
fi
#
gcc -o gmsh_io_prb gmsh_io_prb.o /$HOME/libc/$ARCH/gmsh_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gmsh_io_prb.o."
  exit
fi
#
rm gmsh_io_prb.o
#
./gmsh_io_prb > gmsh_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running gmsh_io_prb."
  exit
fi
rm gmsh_io_prb
#
echo "Program output written to gmsh_io_prb_output.txt"
