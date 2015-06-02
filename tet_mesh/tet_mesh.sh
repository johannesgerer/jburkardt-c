#!/bin/bash
#
cp tet_mesh.h /$HOME/include
#
gcc -c -I /$HOME/include tet_mesh.c
if [ $? -ne 0 ]; then
  echo "Errors compiling tet_mesh.c"
  exit
fi
#
mv tet_mesh.o ~/libc/$ARCH/tet_mesh.o
#
echo "Library installed as ~/libc/$ARCH/tet_mesh.o"
