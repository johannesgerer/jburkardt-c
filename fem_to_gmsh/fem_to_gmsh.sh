#!/bin/bash
#
gcc -c fem_to_gmsh.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_to_gmsh.c"
  exit
fi
#
gcc fem_to_gmsh.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_to_gmsh.o."
  exit
fi
#
rm fem_to_gmsh.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/fem_to_gmsh
#
echo "Executable installed as ~/binc/$ARCH/fem_to_gmsh"
