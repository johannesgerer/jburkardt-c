#!/bin/bash
#
gcc -c flood_opengl.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling flood_opengl.c"
  exit
fi
rm compiler.txt
#
gcc flood_opengl.o -framework OpenGL -framework GLUT
if [ $? -ne 0 ]; then
  echo "Errors linking flood_opengl.o"
  exit
fi
#
rm flood_opengl.o
mv a.out ~/binc/$ARCH/flood_opengl
#
echo "Executable installed as ~/binc/$ARCH/flood_opengl"
