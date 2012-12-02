#!/bin/bash
#
cp mesh_io.h /$HOME/include
#
gcc -c -g -I /$HOME/include mesh_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mesh_io.c."
  exit
fi
rm compiler.txt
#
mv mesh_io.o ~/libc/$ARCH/mesh_io.o
#
echo "Library installed as ~/libc/$ARCH/mesh_io.o"
