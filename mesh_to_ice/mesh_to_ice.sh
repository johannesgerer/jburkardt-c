#!/bin/bash
#
gcc -c -g -I /usr/local/include mesh_to_ice.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling mesh_to_ice.c"
  exit
fi
rm compiler.txt
#
gcc mesh_to_ice.o -L/usr/local/lib -lnetcdf
if [ $? -ne 0 ]; then
  echo "Errors while loading mesh_to_ice.o"
  exit
fi
rm mesh_to_ice.o
#
mv a.out ~/binc/$ARCH/mesh_to_ice
#
echo "Executable installed as ~/binc/$ARCH/mesh_to_ice"
