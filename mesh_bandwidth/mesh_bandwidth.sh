#!/bin/bash
#
gcc -c -g mesh_bandwidth.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mesh_bandwidth.c"
  exit
fi
rm compiler.txt
#
gcc mesh_bandwidth.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mesh_bandwidth.o"
  exit
fi
rm mesh_bandwidth.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/mesh_bandwidth
#
echo "Executable installed as ~/binc/$ARCH/mesh_bandwidth"
