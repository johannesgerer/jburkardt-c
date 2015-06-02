#!/bin/bash
#
cp polygon_triangulate.h /$HOME/include
#
gcc -c -I/$HOME/include polygon_triangulate.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_triangulate.c"
  exit
fi
rm compiler.txt
#
mv polygon_triangulate.o ~/libc/$ARCH/polygon_triangulate.o
#
echo "Library installed as ~/libc/$ARCH/polygon_triangulate.o"
