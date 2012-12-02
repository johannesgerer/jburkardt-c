#!/bin/bash
#
cp mpi_stubs_c.h /$HOME/include/mpi_stubs_c.h
#
gcc -c -g mpi_stubs.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mpi_stubs.c."
  exit
fi
rm compiler.txt
#
mv mpi_stubs.o ~/libc/$ARCH/mpi_stubs.o
#
echo "Library installed as ~/libc/$ARCH/mpi_stubs.o"
