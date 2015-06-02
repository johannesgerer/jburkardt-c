#!/bin/bash
#
cp gmsh_io.h /$HOME/include
#
gcc -c -I/$HOME/include gmsh_io.c
if [ $? -ne 0 ]; then
  echo "Errors compiling gmsh_io.c"
  exit
fi
#
mv gmsh_io.o ~/libc/$ARCH/gmsh_io.o
#
echo "Library installed as ~/libc/$ARCH/gmsh_io.o"
