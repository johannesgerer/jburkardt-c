#!/bin/bash
#
gcc -c ice_to_mesh.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ice_to_mesh.c."
  exit
fi
rm compiler.txt
#
gcc ice_to_mesh.o -L$HOME/lib/$ARCH -lnetcdf -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ice_to_mesh.o."
  exit
fi
#
rm ice_to_mesh.o
#
chmod u+x a.out
mv a.out ~/binc/$ARCH/ice_to_mesh
#
echo "Executable installed as ~/binc/$ARCH/ice_to_mesh"
