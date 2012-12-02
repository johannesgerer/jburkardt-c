#!/bin/bash
#
gcc -c life_opengl.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling life_opengl.c"
  exit
fi
rm compiler.txt
#
gcc life_opengl.o -framework OpenGL -framework GLUT
if [ $? -ne 0 ]; then
  echo "Errors linking life_opengl.o"
  exit
fi
#
rm life_opengl.o
mv a.out ~/binc/$ARCH/life_opengl
#
echo "Executable installed as ~/binc/$ARCH/life_opengl"
