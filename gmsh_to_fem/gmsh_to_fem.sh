#!/bin/bash
#
gcc -c gmsh_to_fem.c
if [ $? -ne 0 ]; then
  echo "Errors compiling gmsh_to_fem.c"
  exit
fi
#
gcc gmsh_to_fem.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gmsh_to_fem.o."
  exit
fi
#
rm gmsh_to_fem.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/gmsh_to_fem
#
echo "Executable installed as ~/binc/$ARCH/gmsh_to_fem"
